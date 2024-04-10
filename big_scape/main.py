"""Main file for BiG-SCAPE cluster module"""

# from python
import sys
import os
import psutil
import signal
import logging
from datetime import datetime
import platform
from pathlib import Path

# from dependencies
import tqdm

# from other modules
from big_scape.cli.config import BigscapeConfig
from big_scape.data import DB
from big_scape.hmm import HMMer, HSP
from big_scape.diagnostics import Profiler
from big_scape.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
    write_record_annotations_file,
    write_full_network_file,
)


import big_scape.file_input as bs_files

# import big_scape.genbank as bs_gbk
import big_scape.data as bs_data
import big_scape.enums as bs_enums
import big_scape.comparison as bs_comparison
import big_scape.network.network as bs_network
import big_scape.network.families as bs_families
import big_scape.distances.mix as bs_mix
import big_scape.distances.legacy_classify as bs_legacy_classify
import big_scape.distances.classify as bs_classify
import big_scape.distances.query as bs_query


def run_bigscape(run: dict) -> None:
    """Run a bigscape cluster analysis. This is the main function of the program that parses the
    command line arguments, loads the data, runs the analysis and saves the output.
    """
    # starting information
    logging.info(
        "Starting BiG-SCAPE %s run on %s level with %s alignment and %s weights",
        run["mode"],
        run["record_type"].value,
        run["alignment_mode"].value,
        "legacy" if run["legacy_weights"] else "mix",
    )

    # root directory
    bigscape_dir = Path(os.path.dirname(os.path.abspath(__file__)))

    main_pid = os.getpid()

    # parse config file
    logging.info("Using config file %s", run["config_file_path"])
    BigscapeConfig.parse_config(run)

    # capture ctrl-c for database saving and other cleanup
    def signal_handler(sig, frame):
        nonlocal main_pid

        # return if not main process
        if os.getpid() != main_pid:
            return
        if DB.opened():
            logging.warning("User requested SIGINT. Saving database and exiting")
            logging.warning("Press ctrl+c again to force exit")

            # kill all child processes
            for child in psutil.Process().children(recursive=True):
                child.kill()

            # reset to default handler in case we get another SIGINT
            signal.signal(signal.SIGINT, signal.SIG_DFL)

            DB.commit()
            DB.save_to_disk(run["db_path"])

        sys.exit(0)

    signal.signal(signal.SIGINT, signal_handler)

    start_time = run["start_time"]

    # at this point a number of things should be non-None. we need to convince mypy of
    # this as well
    assert run["pfam_path"] is not None

    if not HMMer.are_profiles_pressed(run["pfam_path"]):
        logging.warning("HMM files were not pressed!")
        HMMer.press(run["pfam_path"])

    # start profiler
    # TODO: cant run if in macbook, fix to catch this
    if run["profiling"]:
        profiler = Profiler(run["profile_path"])
        profiler.start()
        logging.info("Profiler started")

    # INPUT - get mibig if needed
    if run["mibig_version"]:
        mibig_version_dir = bs_files.get_mibig(run["mibig_version"], bigscape_dir)
        run["mibig_dir"] = mibig_version_dir

    # INPUT - create an in memory DB or load from disk
    if not run["db_path"].exists():
        bs_data.DB.create_in_mem()
    else:
        bs_data.DB.load_from_disk(run["db_path"])

    # INPUT - load data
    gbks = bs_files.load_gbks(run, bigscape_dir)

    # INPUT - get all working BGC records
    if run["query_bgc_path"]:
        all_bgc_records, query_record = bs_files.get_all_bgc_records_query(run, gbks)
    else:
        all_bgc_records = bs_files.get_all_bgc_records(run, gbks)

    # get fist task
    run_state = bs_data.find_minimum_task(gbks)

    # nothing to do!
    if run_state == bs_enums.TASK.NOTHING_TO_DO:
        logging.info("Nothing to do!")
        sys.exit(0)

    logging.info("First task: %s", run_state)

    # we will need this later. this can be set in hmmscan, hmmalign or not at all
    pfam_info = None

    # HMMER - Search/scan
    if run_state == bs_enums.TASK.HMM_SCAN:
        HMMer.init(run["pfam_path"])

        # first opportunity to set this is here
        pfam_info = HMMer.get_pfam_info()

        # we certainly need to do some sort of scan, since the run_state is on hmm_scan. but we do
        # not need to do everything. Find the CDS that need scanning
        cds_to_scan = bs_data.get_cds_to_scan(gbks)

        logging.info("Scanning %d CDS", len(cds_to_scan))

        if platform.system() == "Darwin":
            logging.debug(
                "Running on %s: hmmsearch_simple with %d cores",
                platform.system(),
                run["cores"],
            )
            with tqdm.tqdm(unit="CDS", total=len(HMMer.profiles), desc="HMMSCAN") as t:
                HMMer.hmmsearch_simple(
                    cds_to_scan,
                    domain_overlap_cutoff=run["domain_overlap_cutoff"],
                    cores=run["cores"],
                    callback=lambda x, y: t.update(),
                )
        else:
            logging.debug(
                "Running on %s: hmmsearch_multiprocess with %d cores",
                platform.system(),
                run["cores"],
            )

            with tqdm.tqdm(unit="CDS", total=len(cds_to_scan), desc="HMMSCAN") as t:

                def callback(tasks_done):
                    t.update(tasks_done)

                # TODO: the overlap filtering in this function does not seem to work
                HMMer.hmmsearch_multiprocess(
                    cds_to_scan,
                    domain_overlap_cutoff=run["domain_overlap_cutoff"],
                    cores=run["cores"],
                    callback=callback,
                )

        # TODO: move, or remove after the add_hsp_overlap function is fixed (if it is broken
        # in the first place)
        # this sorts all CDS and then filters them using the old filtering system, which
        # is less efficient than the flitering using the CDS.add_hsp_overlap_filter
        # method. however, that method seems to be broken somehow
        for gbk in gbks:
            for cds in gbk.genes:
                cds.hsps = sorted(cds.hsps)

        exec_time = datetime.now() - start_time
        logging.info("scan done at %f seconds", exec_time.total_seconds())

        # save hsps to database

        for cds in cds_to_scan:
            for hsp in cds.hsps:
                hsp.save(False)
        DB.commit()

        exec_time = datetime.now() - start_time
        logging.info("DB: HSP save done at %f seconds", exec_time.total_seconds())

        HMMer.unload()

        # set new run state
        run_state = bs_data.find_minimum_task(gbks)
        logging.info("Next task: %s", run_state)

        all_hsps: list[HSP] = []
        for gbk in gbks:
            for cds in gbk.genes:
                all_hsps.extend(cds.hsps)

        logging.info("%d hsps found in this run", len(all_hsps))

    # HMMER - Align
    if run_state == bs_enums.TASK.HMM_ALIGN:
        HMMer.init(run["pfam_path"], False)

        # if this wasn't set before, set it now
        if pfam_info is None:
            pfam_info = HMMer.get_pfam_info()

        hsps_to_align = bs_data.get_hsp_to_align(gbks)

        with tqdm.tqdm(unit="HSP", total=len(hsps_to_align), desc="HMMALIGN") as t:

            def align_callback(tasks_done: int):
                t.update(tasks_done)

            HMMer.align_simple(hsps_to_align, align_callback)

        alignment_count = 0

        for hsp in hsps_to_align:
            if hsp.alignment is None:
                continue
            hsp.alignment.save(False)
            alignment_count += 1

        logging.info("%d alignments", alignment_count)

        exec_time = datetime.now() - start_time
        logging.info("align done at %f seconds", exec_time.total_seconds())

        DB.commit()
        DB.save_to_disk(run["db_path"])

        exec_time = datetime.now() - start_time
        logging.info(
            "DB: HSP alignment save done at %f seconds", exec_time.total_seconds()
        )

        HMMer.unload()

        # set new run state
        run_state = bs_data.find_minimum_task(gbks)
        logging.info("Next task: %s", run_state)

    # TODO: idea: use sqlite to set distances of 1.0 for all pairs that have no domains
    # in common

    # DISTANCE GENERATION

    # mix

    if not run["no_mix"] and not run["query_bgc_path"]:
        bs_mix.calculate_distances_mix(run, all_bgc_records)

        DB.commit()

    # legacy_classify

    if run["legacy_classify"] and not run["query_bgc_path"]:
        bs_legacy_classify.calculate_distances_legacy_classify(run, all_bgc_records)

        DB.commit()

    # classify

    if run["classify"] and not run["query_bgc_path"]:
        bs_classify.calculate_distances_classify(run, all_bgc_records)

        DB.commit()

    # query

    if run["query_bgc_path"]:
        query_records = bs_query.calculate_distances_query(
            run, all_bgc_records, query_record
        )

        DB.commit()

    # FAMILY GENERATION
    # TODO: implement include nodes (all bgc records of this run) in cc component
    # & egde param id (account for weights/alignment mode) to use only correct edges
    logging.info("Generating families")

    bs_families.reset_db_families()

    bs_network.reset_db_connected_components()

    # mix

    if not run["no_mix"] and not run["query_bgc_path"]:
        edge_param_id = bs_comparison.get_edge_param_id(run, "mix")
        for cutoff in run["gcf_cutoffs"]:
            logging.debug("Bin 'mix': cutoff %s", cutoff)
            for connected_component in bs_network.get_connected_components(
                cutoff, edge_param_id, all_bgc_records
            ):
                logging.debug(
                    "Found connected component with %d edges", len(connected_component)
                )

                regions_families = bs_families.generate_families(
                    connected_component, "mix", cutoff
                )

                # save families to database
                bs_families.save_to_db(regions_families)
            bs_families.save_singletons(run["record_type"], cutoff, "mix")

        DB.commit()

    # legacy_classify

    if run["legacy_classify"] and not run["query_bgc_path"]:
        for bin in bs_comparison.legacy_bin_generator(all_bgc_records, run):
            edge_param_id = bs_comparison.get_edge_param_id(run, bin.weights)
            for cutoff in run["gcf_cutoffs"]:
                logging.debug("Bin '%s': cutoff %s", bin.label, cutoff)
                for connected_component in bs_network.get_connected_components(
                    cutoff, edge_param_id, bin.source_records
                ):
                    regions_families = bs_families.generate_families(
                        connected_component, bin.label, cutoff
                    )
                    bs_families.save_to_db(regions_families)
                bs_families.save_singletons(run["record_type"], cutoff, bin.label)

        DB.commit()

    # classify

    if run["classify"] and not run["query_bgc_path"]:
        for bin in bs_comparison.as_class_bin_generator(all_bgc_records, run):
            edge_param_id = bs_comparison.get_edge_param_id(run, bin.weights)
            for cutoff in run["gcf_cutoffs"]:
                logging.debug("Bin '%s': cutoff %s", bin.label, cutoff)
                for connected_component in bs_network.get_connected_components(
                    cutoff, edge_param_id, bin.source_records
                ):
                    regions_families = bs_families.generate_families(
                        connected_component, bin.label, cutoff
                    )
                    bs_families.save_to_db(regions_families)
                bs_families.save_singletons(run["record_type"], cutoff, bin.label)

        DB.commit()

    # query BGC mode
    if run["query_bgc_path"]:
        cc_cutoff: dict[
            str, list[tuple[int, int, float, float, float, float, int]]
        ] = {}
        if run["legacy_weights"]:
            weights = bs_comparison.get_legacy_weights_from_category(
                query_record, query_record.product, run
            )
        else:
            weights = "mix"
        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        for cutoff in run["gcf_cutoffs"]:
            logging.info("Query BGC Bin: cutoff %s", cutoff)

            # get_connected_components returns a list of connected components, but we only
            # want the first one, so we use next()

            query_connected_component = next(
                bs_network.get_connected_components(
                    cutoff, edge_param_id, [query_record]
                )
            )

            cc_cutoff[cutoff] = query_connected_component

            logging.debug(
                "Found connected component with %d edges",
                len(query_connected_component),
            )

            regions_families = bs_families.generate_families(
                query_connected_component, weights, cutoff
            )

            # save families to database
            bs_families.save_to_db(regions_families)

        DB.commit()

    DB.save_to_disk(run["db_path"])

    # OUTPUT GENERATION

    exec_time = datetime.now() - start_time
    run["duration"] = exec_time
    run["end_time"] = datetime.now()

    # if this wasn't set in scan or align, set it now
    if pfam_info is None:
        HMMer.init(run["pfam_path"], False)
        pfam_info = HMMer.get_pfam_info()
        HMMer.unload()

    logging.info("Generating GCF alignments, trees and outputfiles")

    # prepare output files
    legacy_prepare_output(run["output_dir"], pfam_info)

    # write full network file
    write_full_network_file(run, all_bgc_records)

    # prepare output files per cutoff
    for cutoff in run["gcf_cutoffs"]:
        # TODO: update to use records and not gbk regions?
        legacy_prepare_cutoff_output(run, cutoff, gbks)
        # write annotations file
        write_record_annotations_file(run, cutoff, all_bgc_records)

    # mix

    if not run["no_mix"] and not run["query_bgc_path"]:
        edge_param_id = bs_comparison.get_edge_param_id(run, "mix")
        mix_bin = bs_comparison.generate_mix_bin(
            all_bgc_records, edge_param_id, run["record_type"]
        )

        for cutoff in run["gcf_cutoffs"]:
            if not run["include_singletons"]:
                mix_bin.cull_singletons(cutoff)
                if len(mix_bin.record_ids) == 0:
                    logging.info(
                        f"Network {mix_bin.label} with cutoff {cutoff} is empty after culling singletons"
                    )
                    continue
            legacy_prepare_bin_output(run, cutoff, mix_bin)
            legacy_generate_bin_output(run, cutoff, mix_bin)

    # legacy_classify

    if run["legacy_classify"] and not run["query_bgc_path"]:
        legacy_class_bins = bs_comparison.legacy_bin_generator(all_bgc_records, run)

        for bin in legacy_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                if not run["include_singletons"]:
                    bin.cull_singletons(cutoff)
                    if len(bin.record_ids) == 0:
                        logging.info(
                            f"Network '{bin.label}' with cutoff {cutoff} is empty after culling singletons"
                        )
                        continue
                legacy_prepare_bin_output(run, cutoff, bin)
                legacy_generate_bin_output(run, cutoff, bin)

    # classify

    if run["classify"] and not run["query_bgc_path"]:
        as_class_bins = bs_comparison.as_class_bin_generator(all_bgc_records, run)

        for bin in as_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                if not run["include_singletons"]:
                    bin.cull_singletons(cutoff)
                    if len(bin.record_ids) == 0:
                        logging.info(
                            f"Network '{bin.label}' with cutoff {cutoff} is empty after culling singletons"
                        )
                        continue
                legacy_prepare_bin_output(run, cutoff, bin)
                legacy_generate_bin_output(run, cutoff, bin)

    # query

    if run["query_bgc_path"]:
        # TODO: select largest cc once, instead of rebuilding bins?
        for cutoff in run["gcf_cutoffs"]:
            query_bin = bs_comparison.ConnectedComponentPairGenerator(
                cc_cutoff[cutoff], "Query"
            )
            query_bin.add_records(
                [record for record in query_records if record is not None]
            )
            legacy_prepare_bin_output(run, cutoff, query_bin)
            legacy_generate_bin_output(run, cutoff, query_bin)

    if run["profiling"]:
        profiler.stop()

    exec_time = datetime.now() - start_time
    run["duration"] = exec_time
    run["end_time"] = datetime.now()

    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
