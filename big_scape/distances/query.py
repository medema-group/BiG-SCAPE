"""contains code for the generation of distances between queries and references in a
given dataset
"""

# TODO: duplicate code with a lot of other edge generation and within this file


# from python
import logging

import tqdm

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison
import big_scape.enums as bs_enums
import big_scape.network.network as bs_network
import big_scape.data as bs_data


def calculate_distances_query(
    run: dict, all_bgc_records: list[bs_gbk.BGCRecord], query_record: bs_gbk.BGCRecord
) -> bs_comparison.RecordPairGenerator:
    """calculates distances between all queries and references in a given dataset and
    saves them to the database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """
    logging.info("Generating query BGC mode bin")

    query_records = get_query_records(run, all_bgc_records, query_record)

    # if legacy weights are on, then use the legacy weights and pass as label to bin generator
    if run["legacy_weights"]:
        weights = bs_comparison.get_legacy_weights_from_category(
            query_record, query_record.product, run
        )
    else:
        weights = "mix"

    max_cutoff = max(run["gcf_cutoffs"])
    edge_param_id = bs_comparison.get_edge_param_id(run, weights)

    query_bin = bs_comparison.QueryRecordPairGenerator(
        weights, edge_param_id, weights, run["record_type"]
    )
    query_bin.add_records(query_records)

    missing_query_bin = bs_comparison.QueryMissingRecordPairGenerator(query_bin)

    calculate_distances(run, missing_query_bin)

    # add last edges

    query_connected_component = next(
        bs_network.get_connected_components(
            max_cutoff, edge_param_id, query_bin, run["run_id"], query_record
        ),
        None,
    )
    if query_connected_component is None:
        # no nodes are connected even with the highest cutoffs in the run
        return query_bin

    query_nodes = bs_network.get_nodes_from_cc(query_connected_component, query_records)

    bs_network.remove_connected_component(
        query_connected_component, query_bin.label, max_cutoff, run["run_id"]
    )

    query_bin_connected = bs_comparison.RecordPairGenerator(
        weights, edge_param_id, weights, record_type=run["record_type"]
    )
    query_bin_connected.add_records(query_nodes)

    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(query_bin_connected)

    calculate_distances(run, missing_edge_bin)

    return query_bin_connected


def get_query_records(
    run: dict, all_bgc_records: list[bs_gbk.BGCRecord], query_record: bs_gbk.BGCRecord
) -> list[bs_gbk.BGCRecord]:
    """returns the query records and checks if the query is a singleton

    Args:
        run (dict): run parameters
        all_bgc_records (list[bs_gbk.BGCRecord]): all records
        query_record (bs_gbk.BGCRecord): the query record

    Returns:
        list[bs_gbk.BGCRecord]: list of query records

    Raises:
        ValueError: if the query record is a singleton
    """

    query_singleton = True
    query_records = [query_record]

    for record in all_bgc_records:
        if record is None:
            continue
        if record.product is None:
            continue
        if record.parent_gbk is None:
            continue
        if record.parent_gbk.source_type == bs_enums.SOURCE_TYPE.QUERY:
            continue

        # if classification mode is off, then use all records
        if not run["classify"]:
            query_singleton = False
            query_records.append(record)

        # if classification mode is on, then use only those records which have the query class/category
        if run["classify"]:
            classify_mode = run["classify"]
            # now we go to all the records, of a given type, of the gbk and check
            # if they have the same class/category as the query
            if classify_mode == bs_enums.CLASSIFY_MODE.CLASS:
                query_class = [query_record.product]
                record_class = [record.product]

                intersect_class = list(set(query_class) & set(record_class))
                if len(intersect_class) > 0:
                    query_singleton = False
                    query_records.append(record)

            if classify_mode == bs_enums.CLASSIFY_MODE.CATEGORY:
                query_category = [query_record.get_category()]
                record_category = [record.get_category()]

                intersect_cats = list(set(query_category) & set(record_category))
                if len(intersect_cats) > 0:
                    query_singleton = False
                    query_records.append(record)

    if query_singleton:
        logging.error(
            f"Query record {query_record.parent_gbk} is a singleton, no other input "
            "records have the same class/category"
        )
        raise ValueError(
            f"Query record {query_record.parent_gbk} is a singleton, no other input "
            "records have the same class/category"
        )

    return query_records


def calculate_distances(run: dict, bin: bs_comparison.RecordPairGenerator):
    """calculates distances between all records in a given dataset and saves them to the
    database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """

    # we will continue propagating until there are no more edges to generate
    while True:
        # fetches the current number of singleton ref <-> connected ref pairs from the database
        num_pairs = bin.num_pairs()

        if num_pairs > 0:
            logging.info("Calculating distances for %d pairs", num_pairs)

            save_batch = []
            num_edges = 0

            with tqdm.tqdm(
                total=num_pairs, unit="edge", desc="Calculating distances"
            ) as t:

                def callback(edges):
                    nonlocal num_edges
                    nonlocal save_batch
                    batch_size = run["cores"] * 100000
                    for edge in edges:
                        num_edges += 1
                        t.update(1)
                        save_batch.append(edge)
                        if len(save_batch) > batch_size:
                            bs_comparison.save_edges_to_db(save_batch, commit=True)
                            save_batch = []

                bs_comparison.generate_edges(
                    bin,
                    run["alignment_mode"],
                    run["extend_strategy"],
                    run["config_file_path"],
                    run["cores"],
                    run["cores"] * 2,
                    callback,
                )

            bs_comparison.save_edges_to_db(save_batch, commit=True)

            bs_data.DB.commit()

            logging.info("Generated %d edges", num_edges)

        if not run["propagate"]:
            # in this case we only want one iteration, the Query -> Ref edges
            break

        if isinstance(bin, bs_comparison.MissingRecordPairGenerator):
            # in this case we only need edges within one cc, no cycles needed
            break

        if isinstance(bin, bs_comparison.QueryMissingRecordPairGenerator):
            # use the num_pairs from the parent bin because in a partial database,
            # all distances for the first cycle(s) might already be present.
            # we still only want to stop when no other connected nodes are discovered.
            if bin.bin.num_pairs() == 0:
                break

            bin.cycle_records(max(run["gcf_cutoffs"]))
