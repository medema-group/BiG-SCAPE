"""contains code for the generation of distances between queries and references in a
given dataset
"""

# from python
import logging

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison
import big_scape.enums as bs_enums
import big_scape.network.network as bs_network


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

        # if classification mode is on, then us only those records which have the query class/category
        if run["classify"]:
            classify_mode = run["classify"]
            # now we go to all the records, of a given type, of the gbk and check
            # if they have the same class/category as the query
            if classify_mode == bs_enums.CLASSIFY_MODE.CLASS:
                query_class = [query_record.product]
                record_class = [record.product]

                if run["hybrids_off"]:
                    query_class = query_class[0].split(".")
                    record_class = record_class[0].split(".")

                intersect_class = list(set(query_class) & set(record_class))
                if len(intersect_class) > 0:
                    query_singleton = False
                    query_records.append(record)

            if classify_mode == bs_enums.CLASSIFY_MODE.CATEGORY:
                query_category = [bs_comparison.get_record_category(query_record)]
                record_category = [bs_comparison.get_record_category(record)]

                if run["hybrids_off"]:
                    query_category = query_category[0].split(".")
                    record_category = record_category[0].split(".")

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
        weights = bs_comparison.get_weight_category(query_record)
    else:
        weights = "mix"

    edge_param_id = bs_comparison.get_edge_param_id(run, weights)

    # generate inital query -> ref pairs
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

        query_edges = bs_comparison.generate_edges(
            missing_edge_bin,
            run["alignment_mode"],
            run["cores"],
            run["cores"] * 2,
        )

        num_edges = 0

        for edge in query_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge)

        logging.info("Generated %d edges", num_edges)

    if run["skip_propagation"]:
        return query_records

    # now we expand these edges from reference to other reference
    # TODO: see if we can implement missing for these
    ref_to_ref_bin = bs_comparison.RefToRefRecordPairGenerator(
        "Ref_Ref", edge_param_id, weights
    )
    ref_to_ref_bin.add_records(query_records)

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

        def callback(done_pairs):
            if num_pairs > 10:
                mod = round(num_pairs / 10)
            else:
                mod = 1

            if done_pairs % mod == 0:
                logging.info(
                    "%d/%d (%.2f%%)",
                    done_pairs,
                    num_pairs,
                    done_pairs / num_pairs * 100,
                )

        # generate edges for these pairs
        ref_edges = bs_comparison.generate_edges(
            ref_to_ref_bin,
            run["alignment_mode"],
            run["cores"],
            run["cores"] * 2,
            callback,
        )

        num_edges = 0
        for edge in ref_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge, True)
        # if no new edges were generated, then break and exit
        # should never happen since generate edges calls num_pairs() which
        # would have broken the loop earlier on.
        if num_edges == 0:
            break

        logging.info("Generated %d edges", num_edges)

    # now we make any last connected ref <-> connected ref pairs that are missing
    # get all the edges in the query connected component
    query_connected_component = bs_network.get_query_connected_component(
        query_records, query_record._db_id, edge_param_id, 1
    )

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

        query_edges = bs_comparison.generate_edges(
            missing_ref_edge_bin, run["alignment_mode"], run["cores"], run["cores"] * 2
        )

        num_edges = 0

        for edge in query_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge)

        logging.info("Generated %d edges", num_edges)

    return query_records
