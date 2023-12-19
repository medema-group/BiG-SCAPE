"""contains code for the generation of distances between all records in a given dataset"""

# from python
import logging

# from dependencies
import tqdm

# from other modules
import big_scape.data as bs_data
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison


def calculate_distances_mix(
    run: dict, list_bgc_records: list[bs_gbk.BGCRecord]
) -> None:
    """calculates distances between all records in a given dataset and saves them to the
    database

    Args:
        run (bs_param.RunParameters): run parameters
        gbks (list[bs_gbk.GBK]): list of gbks
    """

    logging.info("Generating mix bin")

    edge_param_id = bs_comparison.get_edge_param_id(run, "mix")

    mix_bin = bs_comparison.generate_mix_bin(list_bgc_records, edge_param_id)

    logging.info(mix_bin)

    # check if there are existing distances
    missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)
    num_pairs = missing_edge_bin.num_pairs()

    if num_pairs > 0:
        logging.info("Calculating distances for %d pairs", num_pairs)

        save_batch = []
        num_edges = 0

        with tqdm.tqdm(total=num_pairs, unit="edge", desc="Calculating distances") as t:

            def callback(edges):
                nonlocal num_edges
                nonlocal save_batch
                batch_size = run["cores"] * 100000
                for edge in edges:
                    num_edges += 1
                    t.update(1)
                    save_batch.append(edge)
                    if len(save_batch) > batch_size:
                        bs_comparison.save_edges_to_db(save_batch)
                        bs_data.DB.commit()
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
