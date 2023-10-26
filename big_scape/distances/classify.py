"""contains code for the generation of antiSMASH class bins and calculating
 distances between all records in a given bin"""

# from python
import logging

# from dependencies
import tqdm

# from other modules
import big_scape.data as bs_data
import big_scape.genbank as bs_gbk
import big_scape.comparison as bs_comparison
import big_scape.enums as bs_enums


def calculate_distances_classify(
    run: dict, gbks: list[bs_gbk.GBK], classify_mode: bs_enums.CLASSIFY_MODE
) -> None:
    """
    calculates distances between all records in a given antismash class bin, based on either
     mix weights or legacy weights and saves them to the database
    """

    logging.info("Generating antismash class bins")

    if run["legacy_weights"]:
        weight_type = "legacy_weights"
    else:
        weight_type = "mix"

    all_records: list[bs_gbk.BGCRecord] = []

    for gbk in gbks:
        if gbk.region is None:
            continue
        if gbk.region.product is None:
            continue
        gbk_records = bs_gbk.bgc_record.get_sub_records(gbk.region, run["record_type"])
        all_records.extend(gbk_records)

    as_class_bins = bs_comparison.as_class_bin_generator(
        all_records, weight_type, classify_mode
    )

    for bin in as_class_bins:
        logging.info(bin)

        # check if there are existing distances
        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(bin)
        num_pairs = missing_edge_bin.num_pairs()

        if num_pairs > 0:
            logging.info("Calculating distances for %d pairs", num_pairs)

            last = 0

            def callback(done_pairs):
                nonlocal last

                if done_pairs - last > 10000:
                    last = done_pairs

            as_class_edges = bs_comparison.generate_edges(
                missing_edge_bin, run["alignment_mode"], run["cores"], callback
            )

            with tqdm.tqdm(
                total=num_pairs, unit="edge", desc="Calculating distances"
            ) as t:
                num_edges = 0
                save_batch = []
                batch_size = run["cores"] * 100000
                for edge in as_class_edges:
                    num_edges += 1
                    t.update(1)
                    save_batch.append(edge)
                    if len(save_batch) > batch_size:
                        bs_comparison.save_edges_to_db(save_batch)
                        bs_data.DB.commit()
                        save_batch = []

            bs_comparison.save_edges_to_db(save_batch)

            bs_data.DB.commit()

            logging.debug(
                f"Generated {num_edges} edges for {bin.label} with weights {bin.weights}"
            )
