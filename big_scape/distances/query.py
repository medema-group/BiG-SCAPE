"""contains code for the generation of distances between queries and references in a
given dataset
"""

# from python
import logging

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison

# import big_scape.enums as bs_enums


def calculate_distances_query(run: dict, gbks: list[bs_gbk.GBK]) -> None:
    """calculates distances between all queries and references in a given dataset and
    saves them to the database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """
    logging.info("Generating query BGC mode bin")

    query_bgcs_records: list[bs_gbk.BGCRecord] = []

    # get all working records
    for gbk in gbks:
        if gbk.region is None:
            continue
        if gbk.region.product is None:
            continue
        query_bgcs_records.append(gbk.region)

    #     # currently query only works for regions, since users have no way to
    #     # specify a protocluster/core in a region that has multiple
    #     if gbk.source_type == bs_enums.SOURCE_TYPE.QUERY:
    #         query_record = bs_gbk.bgc_record.get_sub_records(gbk.region, "region")
    #         query_bgcs_records.extend(query_record)
    #     if not run["classify"]:
    #         gbk_records = bs_gbk.bgc_record.get_sub_records(gbk.region, "region")
    #         query_bgcs_records.extend(gbk_records)

    # # if legacy weights are on, then use the legacy weights and pass as label to bin generator
    # if run["legacy_weights"]:
    #     weights = bs_comparison.get_weight_category(query_record)
    # else:
    #     weights = "mix"

    # # if classification mode is on, then us only those records which have the query class/category
    # if run["classify"]:
    #     classify_mode = run["classify"]

    #     if classify_mode == bs_enums.CLASSIFY_MODE.CLASS:
    #         query_class = gbk.region.product
    #         for record in query_bgcs_records:
    #             record_class = record.product
    #             if record_class == query_class:
    #                 query_bgcs_records.extend(record)

    #     if classify_mode == bs_enums.CLASSIFY_MODE.CATEGORY:
    #         query_category = bs_comparison.get_record_category(gbk.region)
    #         for record in query_bgcs_records:
    #             record_category = bs_comparison.get_record_category(record)
    #             if record_category == query_category:
    #                 query_bgcs_records.extend(record)

    # generate inital query -> ref pairs
    query_to_ref_bin = bs_comparison.QueryToRefRecordPairGenerator(
        label="Query_Ref", weights="mix"
    )
    query_to_ref_bin.add_records(query_bgcs_records)

    # fetch any existing distances from database
    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(query_to_ref_bin)
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
        )

        num_edges = 0

        for edge in query_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge)

        logging.info("Generated %d edges", num_edges)

    # now we expand these edges from reference to other reference
    # TODO: see if we can implement missing for these
    ref_to_ref_bin = bs_comparison.RefToRefRecordPairGenerator(
        label="Ref_Ref", weights="mix"
    )
    ref_to_ref_bin.add_records(query_bgcs_records)

    while True:
        # fetches the current number of singleton ref <-> connected ref pairs from the database
        num_pairs = ref_to_ref_bin.num_pairs()
        # if there are no more singleton ref <-> connected ref pairs, then break and exit
        if num_pairs == 0:
            break

        logging.info(
            "Calculating distances for %d pairs (Reference to Reference)",
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
            ref_to_ref_bin, run["alignment_mode"], run["cores"], callback
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
