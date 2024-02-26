"""contains code for the generation of antiSMASH class bins and calculating
 distances between all records in a given bin"""

# TODO: duplicate code with a lot of other edge generation

# from python
import logging

# from dependencies
import tqdm

# from other modules
import big_scape.data as bs_data
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison


def calculate_distances_classify(
    run: dict, all_records: list[bs_gbk.BGCRecord]
) -> None:
    """
    calculates distances between all records in a given antismash class bin, based on either
     mix weights or legacy weights and saves them to the database
    """

    logging.info("Generating antismash class bins")

    as_class_bins = bs_comparison.as_class_bin_generator(all_records, run)

    for bin in as_class_bins:
        logging.info(bin)

        # check if there are existing distances
        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(bin)
        num_pairs = missing_edge_bin.num_pairs()

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
                    missing_edge_bin,
                    run["alignment_mode"],
                    run["cores"],
                    run["cores"] * 2,
                    callback,
                )

            bs_comparison.save_edges_to_db(save_batch)

            bs_data.DB.commit()

            logging.debug(
                f"Generated {num_edges} edges for {bin.label} with weights {bin.weights}"
            )
