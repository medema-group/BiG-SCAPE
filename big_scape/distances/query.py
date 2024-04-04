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
) -> list[bs_gbk.BGCRecord]:
    """calculates distances between all queries and references in a given dataset and
    saves them to the database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """
    logging.info("Generating query BGC mode bin")

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
                query_category = [bs_comparison.get_record_category(query_record)]
                record_category = [bs_comparison.get_record_category(record)]

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

    # if legacy weights are on, then use the legacy weights and pass as label to bin generator
    if run["legacy_weights"]:
        weights = bs_comparison.get_legacy_weights_from_category(
            query_record, query_record.product, run
        )
    else:
        weights = "mix"

    edge_param_id = bs_comparison.get_edge_param_id(run, weights)

    # generate initial query -> ref pairs
    query_to_ref_bin = bs_comparison.QueryToRefRecordPairGenerator(
        "Query_Ref", edge_param_id, weights
    )
    query_to_ref_bin.add_records(query_records)

    # fetch any existing distances from database
    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(query_to_ref_bin)
    # get the number of pairs that are missing
    num_pairs = missing_edge_bin.num_pairs()

    # calculate distances
    if num_pairs > 0:
        logging.info(
            "Calculating distances for %d pairs (Query to Reference)",
            num_pairs,
        )

        save_batch = []
        num_edges = 0

        with tqdm.tqdm(
            total=num_pairs, unit="edge", desc="Calculating distances (qry2ref)"
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
                missing_edge_bin,
                run["alignment_mode"],
                run["cores"],
                run["cores"] * 2,
                callback,
            )

        bs_comparison.save_edges_to_db(save_batch)

        bs_data.DB.commit()

        logging.info("Generated %d edges", num_edges)

    if run["skip_propagation"]:
        return query_records

    # at this point we have query -> ref edges
    # we now need to propagate edges from those ref -> other ref

    # TODO: see if we can implement missing for these
    ref_to_ref_bin = bs_comparison.RefToRefRecordPairGenerator(
        "Ref_Ref", edge_param_id, weights
    )
    ref_to_ref_bin.add_records(query_records)

    # we will continue propagating until there are no more edges to generate
    while True:
        # fetches the current number of singleton ref <-> connected ref pairs from the database
        num_pairs = ref_to_ref_bin.num_pairs()

        # if there are no more singleton ref <-> connected ref pairs, then break and exit
        if num_pairs == 0:
            break

        logging.info(
            "Calculating distances for %d pairs (Connected Reference to Singleton Reference)",
            num_pairs,
        )
        save_batch = []
        num_edges = 0

        with tqdm.tqdm(
            total=num_pairs, unit="edge", desc="Calculating distances (ref2sng)"
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

            # generate edges for these pairs
            bs_comparison.generate_edges(
                ref_to_ref_bin,
                run["alignment_mode"],
                run["cores"],
                run["cores"] * 2,
                callback,
            )

        bs_comparison.save_edges_to_db(save_batch)

        logging.info("Generated %d edges", num_edges)

    logging.info("Finished propagating edges")

    # now we make any last connected ref <-> connected ref pairs that are missing
    # get all the edges in the query connected component
    query_connected_component = bs_network.get_connected_components(
        1, edge_param_id, [query_record]
    )

    # get_connected_components returns a list of connected components, we only want the first one
    query_connected_component = next(query_connected_component)

    query_nodes = bs_network.get_nodes_from_cc(query_connected_component, query_records)

    query_connected_bin = bs_comparison.RecordPairGenerator(
        "Query", edge_param_id, record_type=run["record_type"]
    )
    query_connected_bin.add_records(query_nodes)

    # fetch any existing distances from database
    missing_ref_edge_bin = bs_comparison.MissingRecordPairGenerator(query_connected_bin)
    # get the number of pairs that are missing
    num_pairs = missing_ref_edge_bin.num_pairs()

    # calculate distances
    if num_pairs > 0:
        logging.info(
            "Calculating distances for %d pairs (Connected Reference to Connected Reference)",
            num_pairs,
        )

        save_batch = []
        num_edges = 0

        with tqdm.tqdm(
            total=num_pairs, unit="edge", desc="Calculating distances (ref2ref)"
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
                missing_ref_edge_bin,
                run["alignment_mode"],
                run["cores"],
                run["cores"] * 2,
                callback,
            )

        bs_comparison.save_edges_to_db(save_batch)

        bs_data.DB.commit()

        logging.info("Generated %d edges", num_edges)

    return query_records
