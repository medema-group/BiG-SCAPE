"""contains code for the generation of distances between queries and references in a
given dataset
"""

# from python
import logging

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison


def calculate_distances_query(run: dict, gbks: list[bs_gbk.GBK]) -> None:
    """calculates distances between all queries and references in a given dataset and
    saves them to the database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """
    logging.info("Generating query BGC mode bin")

    query_bgcs_records: list[bs_gbk.BGCRecord] = []

    for gbk in gbks:
        if gbk.region is None:
            continue
        query_bgcs_records.append(gbk.region)

    query_to_ref_bin = bs_comparison.QueryToRefRecordPairGenerator("mix")
    query_to_ref_bin.add_records(query_bgcs_records)

    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(query_to_ref_bin)
    num_pairs = missing_edge_bin.num_pairs()

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
    ref_to_ref_bin = bs_comparison.RefToRefRecordPairGenerator("Ref_Ref")
    ref_to_ref_bin.add_records(query_bgcs_records)

    while True:
        num_pairs = ref_to_ref_bin.num_pairs()

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

        ref_edges = bs_comparison.generate_edges(
            ref_to_ref_bin, run["alignment_mode"], run["cores"], callback
        )

        num_edges = 0
        for edge in ref_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge, True)

        if num_edges == 0:
            break

        logging.info("Generated %d edges", num_edges)
