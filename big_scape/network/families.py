"""Contains functions to execute affinity propagaion on arbitrary graphs"""

# from python
import sys
from typing import Callable, Optional
import warnings
import numpy as np
import networkx
import math
import logging

# from dependencies
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning
from sqlalchemy import select

# from other modules
from big_scape.data import DB
from big_scape.cli.config import BigscapeConfig
from big_scape.genbank.bgc_record import BGCRecord
import big_scape.network.network as bs_network
import big_scape.comparison as bs_comparison

# from this module
from .utility import edge_list_to_sim_matrix


def generate_families(
    connected_component: list[tuple[int, int, float, float, float, float, int]],
    bin_label: str,
    cutoff: float,
    run_id: int,
) -> list[tuple[int, int, float, str, int]]:
    """Execute affinity propagation on a connected component

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, str]]):
            connected component in the form of a list of edges
        bin_label (str): name of current bin to generate families for
        cutoff (float): cutoff used in generation of the connected_component
        run_id (int): id of the current run

    Returns:
        list[tuple[int, int, float, int, int]]: list of
            (region_id, family, cutoff, edge_param_id, run_id) tuples
    """
    # assemble list of (region_id, family, cutoff, edge_param_id, run_id) tuples for easy
    # insertion into db
    regions_families = []

    similarity_matrix, node_ids = edge_list_to_sim_matrix(connected_component)

    if get_cc_density(connected_component) > BigscapeConfig.DENSITY:
        # if a connected component is highly connected, no (or less) splitting is needed
        # run affinity propagation with a lower preference to find the best family center
        labels, centers = aff_sim_matrix(
            similarity_matrix, BigscapeConfig.DENSE_PREFERENCE
        )
    else:
        labels, centers = aff_sim_matrix(similarity_matrix)

    # If affinity propagation did not converge, no centers are returned.
    # to show them in the network anyways, merge them into one arbitrary family
    if len(centers) == 0:
        center = node_ids[0]

        if DB.metadata is None:
            raise RuntimeError("DB metadata is None!")

        gbk_table = DB.metadata.tables["gbk"]
        record_table = DB.metadata.tables["bgc_record"]
        center_data = DB.execute(
            select(
                gbk_table.c.path,
                record_table.c.record_type,
                record_table.c.record_number,
            )
            .join(record_table, record_table.c.gbk_id == gbk_table.c.id)
            .where(record_table.c.id == center)
        ).fetchone()

        if center_data is None:
            raise RuntimeError("Family center not found in database: %s", center)

        c_path, c_type, c_number = center_data
        logging.warning(
            "Affinity Propagation did not converge, records in this connected component "
            "have been merged into one arbitrary family with center: %s",
            "_".join(map(str, [c_path.split("/")[-1], c_type, c_number])),
        )
        return [(rec_id, center, cutoff, bin_label, run_id) for rec_id in node_ids]

    for idx, label in enumerate(labels):
        label = int(label)
        if label == -1:
            continue

        region_id = node_ids[idx]

        center = centers[label]
        family = node_ids[center]

        regions_families.append((region_id, family, cutoff, bin_label, run_id))

    return regions_families


def get_cc_edge_weight_std(connected_component) -> float:
    """calculates the standard deviation of the edge weights of a connected component

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, str]]):
            connected component in the form of a list of edges

    Returns:
        float: standard deviation of the edge weights of the connected component
    """

    edge_weights = [edge[2] for edge in connected_component]
    edge_std = np.std(edge_weights)
    edge_std = round(edge_std, 2)

    return edge_std


def get_cc_density(
    connected_component: list[tuple[int, int, float, float, float, float, int]]
) -> float:
    """calculates the density of a connected component: nr edges / nr of possible edges

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, int]]):
            connected component in the form of a list of edges

    Returns:
        float: density of the connected component
    """

    nr_edges = len(connected_component)

    nodes_a = [edge[0] for edge in connected_component]
    nodes_b = [edge[1] for edge in connected_component]
    nr_nodes = len(set(nodes_a + nodes_b))

    cc_density = nr_edges / (nr_nodes * (nr_nodes - 1) / 2)
    cc_density = round(cc_density, 2)

    return cc_density


def test_centrality(connected_component, node_fraction) -> tuple[bool, list[int]]:
    """tests if a network will break when removing the top nodes
    with highest betweenness centrality

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, str]]):
            connected component in the form of a list of edges
        node_fraction (float): fraction of nodes with highest betweenness centrality to remove

    Returns:
        tuple[bool, list[int]]: whether the network breaks and the list of nodes sorted by betweenness centrality
    """

    edgelist = [(edge[0], edge[1], edge[2]) for edge in connected_component]

    graph = networkx.Graph()
    graph.add_weighted_edges_from(edgelist)

    betweeness_centrality_dict = networkx.betweenness_centrality(graph)
    sorted_between_bentrality_nodes = sorted(
        betweeness_centrality_dict, key=betweeness_centrality_dict.get, reverse=True
    )

    # round up to nearest integer
    top_nodes = math.ceil(len(sorted_between_bentrality_nodes) * node_fraction)
    nodes_to_remove = sorted_between_bentrality_nodes[:top_nodes]

    for node in nodes_to_remove:
        graph.remove_node(node)

    nr_ccs = networkx.number_connected_components(graph)

    del graph

    if nr_ccs > 1:
        return True, sorted_between_bentrality_nodes

    return False, sorted_between_bentrality_nodes


def aff_sim_matrix(matrix, preference: Optional[float] = None):
    """Execute affinity propagation on a __similarity__ matrix

    Note: a similarity matrix. Not a distance matrix.

    Args:
        matrix (numpy.array[numpy.array]): similarity matrix in numpy array of array
        format.
        preference (float, optional): Affinity propagation preference.

    Returns:
        tuple[list[int], list[int]]: list of labels and list of cluster center ids
    """
    # thanks numpy but we sort of know what we're doing
    warnings.filterwarnings(action="ignore", category=ConvergenceWarning)

    if preference is None:
        preference = BigscapeConfig.PREFERENCE

    af_results = AffinityPropagation(
        damping=0.90,
        max_iter=1000,
        convergence_iter=200,
        affinity="precomputed",
        preference=preference,
    ).fit(matrix)

    return af_results.labels_, af_results.cluster_centers_indices_


def save_to_db(regions_families: list[tuple[int, int, float, str, int]]):
    """Save families to database

    Args:
        regions_families (list[tuple[int, int, float, str, int]]): list of (region_id, family,
        cutoff, bin_label, run_id) tuples
    """
    if DB.metadata is None:
        raise RuntimeError("DB metadata is None!")
    family_table = DB.metadata.tables["family"]
    bgc_record_family_table = DB.metadata.tables["bgc_record_family"]

    for region_id, family, cutoff, bin_label, run_id in regions_families:
        # obtain unique family id if present
        fam_id_query = (
            select(family_table.c.id)
            .where(family_table.c.center_id == family)
            .where(family_table.c.cutoff == cutoff)
            .where(family_table.c.bin_label == bin_label)
            .where(family_table.c.run_id == run_id)
        )
        family_id = DB.execute(fam_id_query).fetchone()

        # if not yet present, insert it
        if family_id is None:
            insert_query = (
                family_table.insert()
                .returning(family_table.c.id)
                .values(
                    center_id=family, cutoff=cutoff, bin_label=bin_label, run_id=run_id
                )
                .compile()
            )

            cursor_result = DB.execute(insert_query, False)
            family_id = cursor_result.fetchone()
            if family_id is None:
                raise RuntimeError("No return value from insert query")

        insert_statement = (
            bgc_record_family_table.insert()
            .values(record_id=region_id, family_id=family_id[0])
            .prefix_with("OR REPLACE")
        )

        DB.execute(insert_statement)


def reset_db_family_tables():
    """Clear previous family assignments from database"""
    DB.execute(DB.metadata.tables["bgc_record_family"].delete())
    DB.execute(DB.metadata.tables["family"].delete())


def save_singletons(
    record_ids: list[int], cutoff: float, bin_label: str, run_id: int
) -> None:
    """Create unique family for any singletons and save to database

    Args:
        record_ids (list[int]): record ids to create any singleton families for
        cutoff (float): cutoff value to create families for
        bin_label (str): label of the bin to create families for
        run_id (int): id of the current run
    """
    if DB.metadata is None:
        raise RuntimeError("DB metadata is None!")

    family_table = DB.metadata.tables["family"]
    bgc_record_family_table = DB.metadata.tables["bgc_record_family"]
    record_table = DB.metadata.tables["bgc_record"]

    singleton_query = (
        select(record_table.c.id)
        .where(record_table.c.id.in_(record_ids))
        .where(
            record_table.c.id.not_in(
                select(bgc_record_family_table.c.record_id)
                .join(
                    family_table,
                    family_table.c.id == bgc_record_family_table.c.family_id,
                )
                .where(family_table.c.cutoff == cutoff)
                .where(family_table.c.bin_label == bin_label)
                .where(family_table.c.run_id == run_id)
            )
        )
    )

    singletons = DB.execute(singleton_query)
    if singletons is None:
        return

    singleton_regions = []
    for singleton in singletons.fetchall():
        singleton_regions.append(
            (singleton[0], singleton[0], cutoff, bin_label, run_id)
        )

        DB.execute(
            DB.metadata.tables["connected_component"]
            .insert()
            .values(
                id=singleton[0],
                record_id=singleton[0],
                cutoff=cutoff,
                bin_label=bin_label,
                run_id=run_id,
            )
        )

    save_to_db(singleton_regions)


def run_family_assignments(
    run: dict, bin_generator: Callable, all_bgc_records: list[BGCRecord]
) -> None:
    """Run the family assignment workflow

    Args:
        run (dict): run configuration
        bin_generator (function): generator function for bins
        all_bgc_records (list[BgcRecord]): all BGC records

    Returns:
        None
    """

    bins = bin_generator(all_bgc_records, run)

    # in the case of the mix bin generator, the function returns a single
    # RecordPairGenerator object, in the classify case, it
    # returns an Iterator of RecordPairGenerator objects
    if isinstance(bins, bs_comparison.RecordPairGenerator):
        bins = [bins]

    for bin in bins:
        edge_param_id = bs_comparison.get_edge_param_id(run, bin.weights)

        for cutoff in run["gcf_cutoffs"]:
            logging.info(
                "Generating connected components for Bin '%s': cutoff %s",
                bin.label,
                cutoff,
            )

            for connected_component in bs_network.get_connected_components(
                cutoff, edge_param_id, bin, run["run_id"]
            ):
                # check and remove ref only cc
                if bs_network.reference_only_connected_component(
                    connected_component, bin.source_records
                ):
                    bs_network.remove_connected_component(
                        connected_component, cutoff, run["run_id"]
                    )
                    continue

                logging.debug(
                    "Found connected component with %d edges",
                    len(connected_component),
                )
                regions_families = generate_families(
                    connected_component, bin.label, cutoff, run["run_id"]
                )
                save_to_db(regions_families)

            if run["include_singletons"]:
                save_singletons(
                    bin.get_query_source_record_ids(), cutoff, bin.label, run["run_id"]
                )

    DB.commit()


def run_family_assignments_query(
    run: dict, query_bin: bs_comparison.RecordPairGenerator, query_record: BGCRecord
) -> dict:
    """Run the family assignment workflow for a single query record

    Args:
        run (dict): run configuration
        query_bin (bs_comparison.RecordPairGenerator): query bin, as generated from
        distance calculations, contains all records at cutoff 1
        query_record (BGCRecord): query record

    """

    cc_cutoff: dict[str, list[tuple[int, int, float, float, float, float, int]]] = {}

    for cutoff in run["gcf_cutoffs"]:
        logging.info("Query BGC Bin: cutoff %s", cutoff)

        # get_connected_components returns a list of connected components, but we only
        # want the first one, so we use next()

        query_connected_component = next(
            bs_network.get_connected_components(
                cutoff,
                query_bin.edge_param_id,
                query_bin,
                run["run_id"],
                query_record,
            ),
            None,
        )

        if query_connected_component is None:
            logging.warning(
                "No connected components found for %s bin at cutoff %s",
                query_bin.label,
                cutoff,
            )
            continue

        cc_cutoff[cutoff] = query_connected_component

        logging.debug(
            "Found connected component with %d edges",
            len(query_connected_component),
        )

        regions_families = generate_families(
            query_connected_component, query_bin.label, cutoff, run["run_id"]
        )

        # save families to database
        save_to_db(regions_families)

    DB.commit()

    # no connected components found
    if cc_cutoff == {}:
        logging.warning(
            "No connected components found for %s bin, stopping run. "
            "The edges generated in this run are still saved to the database,"
            " so you can try and re-running with a less strict cutoff(s)",
            query_bin.label,
        )
        sys.exit(0)

    return cc_cutoff
