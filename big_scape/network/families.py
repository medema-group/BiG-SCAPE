"""Contains functions to execute affinity propagaion on arbitrary graphs"""

# from python
import warnings
import numpy as np
import networkx
import math

# from dependencies
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning
from sqlalchemy import select

# from other modules
from big_scape.data import DB
from big_scape.enums import RECORD_TYPE
from big_scape.cli.config import BigscapeConfig

# from this module
from .utility import edge_list_to_sim_matrix


def generate_families(
    connected_component: list[tuple[int, int, float, float, float, float, int]],
    bin_label: str,
    cutoff: float,
) -> list[tuple[int, int, float, str]]:
    """Execute affinity propagation on a connected component

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, str]]):
            connected component in the form of a list of edges
        bin_label (str): name of current bin to generate families for
        cutoff (float): cutoff used in generation of the connected_component

    Returns:
        list[tuple[int, int, float, int]]: list of
            (region_id, family, cutoff, edge_param_id) tuples
    """
    # assemble list of (region_id, family, cutoff, edge_param_id) tuples for easy
    # insertion into db
    regions_families = []

    # hierarchical checks to see if affinity propagation should be applied or not
    # AP is applied if the edge weight standard deviation is above a threshold
    # or if the connected component connectivity is below a threshold
    # or if the connected component breaks when removing the top nodes with
    # highest betweenness centrality
    # TODO: consider getting rid of centrality check, but still
    # need a center node for trees
    cc_edge_weight_std = get_cc_edge_weight_std(connected_component)
    if cc_edge_weight_std < BigscapeConfig.EDGE_WEIGHT_STD_THRESHOLD:
        cc_connectivity = get_cc_connectivity(connected_component)
        if cc_connectivity > BigscapeConfig.CC_CONNECTIVITY_THRESHOLD:
            cc_breaks, nodes_centrality = test_centrality(
                connected_component, BigscapeConfig.BETWEENNESS_CENTRALITY_NODES
            )
            if not cc_breaks:
                # node with highest betweenness centrality is the center of the family
                family_id = nodes_centrality[0]

                for edge in connected_component:
                    regions_families.append((edge[0], family_id, cutoff, bin_label))
                    regions_families.append((edge[1], family_id, cutoff, bin_label))

                return regions_families

    distance_matrix, node_ids = edge_list_to_sim_matrix(connected_component)

    labels, centers = aff_sim_matrix(distance_matrix)

    for idx, label in enumerate(labels):
        label = int(label)
        if label == -1:
            continue

        region_id = node_ids[idx]

        center = centers[label]
        family = node_ids[center]

        regions_families.append((region_id, family, cutoff, bin_label))

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


def get_cc_connectivity(connected_component) -> float:
    """calculates the connectivity of a connected component

    Args:
        connected_component (list[tuple[int, int, float, float, float, float, str]]):
            connected component in the form of a list of edges

    Returns:
        float: connectivity of the connected component
    """

    nr_edges = len(connected_component)

    nodes_a = [edge[0] for edge in connected_component]
    nodes_b = [edge[1] for edge in connected_component]
    nr_nodes = len(set(nodes_a + nodes_b))

    cc_connectivity = nr_edges / (nr_nodes * (nr_nodes - 1) / 2)
    cc_connectivity = round(cc_connectivity, 2)

    return cc_connectivity


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


def aff_sim_matrix(matrix):
    """Execute affinity propagation on a __similarity__ matrix

    Note: a similarity matrix. Not a distance matrix.

    Args:
        matrix (numpy.array[numpy.array]): similarity matrix in numpy array of array
        format.

    Returns:
        tuple[list[int], list[int]]: list of labels and list of cluster center ids
    """
    # thanks numpy but we sort of know what we're doing
    warnings.filterwarnings(action="ignore", category=ConvergenceWarning)

    af_results = AffinityPropagation(
        damping=0.90,
        max_iter=1000,
        convergence_iter=200,
        affinity="precomputed",
    ).fit(matrix)

    return af_results.labels_, af_results.cluster_centers_indices_


def save_to_db(regions_families):
    """Save families to database

    Args:
        regions_families (list[tuple[int, int, float]]): list of (region_id, family,
        cutoff) tuples
    """
    family_table = DB.metadata.tables["family"]
    bgc_record_family_table = DB.metadata.tables["bgc_record_family"]

    for region_id, family, cutoff, bin_label in regions_families:
        # obtain unique family id if present
        fam_id_query = (
            select(family_table.c.id)
            .where(family_table.c.center_id == family)
            .where(family_table.c.cutoff == cutoff)
            .where(family_table.c.bin_label == bin_label)
        )
        family_id = DB.execute(fam_id_query).fetchone()

        # if not yet present, insert it
        if family_id is None:
            insert_query = (
                family_table.insert()
                .returning(family_table.c.id)
                .values(center_id=family, cutoff=cutoff, bin_label=bin_label)
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


def reset_db_families():
    """Clear previous family assignments from database"""
    DB.execute(DB.metadata.tables["bgc_record_family"].delete())
    DB.execute(DB.metadata.tables["family"].delete())


def save_singletons(record_type: RECORD_TYPE, cutoff: float, bin_label: str) -> None:
    """Create unique family for any singletons and save to database

    Args:
        record_type (RECORD_TYPE): record type to create families for
        cutoff (float): cutoff value to create families for
        bin_label (str): label of the bin to create families for
    """
    record_type_str = record_type.value

    if DB.metadata is None:
        raise RuntimeError("DB metadata is None!")

    family_table = DB.metadata.tables["family"]
    bgc_record_family_table = DB.metadata.tables["bgc_record_family"]
    record_table = DB.metadata.tables["bgc_record"]

    singleton_query = (
        select(record_table.c.id)
        .where(record_table.c.record_type == record_type_str)
        .where(
            record_table.c.id.not_in(
                select(bgc_record_family_table.c.record_id)
                .join(
                    family_table,
                    family_table.c.id == bgc_record_family_table.c.family_id,
                )
                .where(family_table.c.cutoff == cutoff)
                .where(family_table.c.bin_label == bin_label)
            )
        )
    )

    singletons = DB.execute(singleton_query)
    if singletons is None:
        return

    singleton_regions = []
    for singleton in singletons.fetchall():
        singleton_regions.append((singleton[0], singleton[0], cutoff, bin_label))

    save_to_db(singleton_regions)
