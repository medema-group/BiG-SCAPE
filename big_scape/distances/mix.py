"""contains code for the generation of distances between all records in a given dataset"""

# from python
import logging

# from other modules
import big_scape.parameters as bs_param
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison


def calculate_distances_mix(
    run: bs_param.RunParameters, gbks: list[bs_gbk.GBK]
) -> None:
    """calculates distances between all records in a given dataset and saves them to the
    database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """

    logging.info("Generating mix bin")

    mix_bgc_regions: list[bs_gbk.BGCRecord] = []

    for gbk in gbks:
        if gbk.region is not None:
            mix_bgc_regions.append(gbk.region)

    mix_bin = bs_comparison.RecordPairGenerator("mix")
    mix_bin.add_records(mix_bgc_regions)

    logging.info(mix_bin)

    # check if there are existing distances
    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)
    num_pairs = missing_edge_bin.num_pairs()

    if num_pairs > 0:
        logging.info("Calculating distances for %d pairs", num_pairs)

        def callback(done_pairs):
            if mix_bin.num_pairs() > 10:
                mod = round(mix_bin.num_pairs() / 10)
            else:
                mod = 1

            if done_pairs % mod == 0:
                logging.info(
                    f"{done_pairs}/{mix_bin.num_pairs()} ({done_pairs/mix_bin.num_pairs():.2%})"
                )

        mix_edges = bs_comparison.generate_edges(
            missing_edge_bin, run.comparison.alignment_mode, run.cores, callback
        )

        num_edges = 0
        for edge in mix_edges:
            num_edges += 1
            bs_comparison.save_edge_to_db(edge)

        logging.info("Generated %d edges", num_edges)
