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
from big_scape.data import DB
from big_scape.hmm import HMMer, HSP
from big_scape.diagnostics import Profiler
from big_scape.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)


import big_scape.file_input as bs_files
import big_scape.genbank as bs_gbk
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
    # root directory
    bigscape_dir = Path(os.path.dirname(os.path.abspath(__file__)))

    main_pid = os.getpid()

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

    # INPUT - load data

    gbks = bs_files.load_gbks(run, bigscape_dir)

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
            logging.warning("Running on mac-OS: hmmsearch_simple single threaded")
            HMMer.hmmsearch_simple(cds_to_scan, 1)
        else:
            logging.debug(
                "Running on %s: hmmsearch_multiprocess with %d cores",
                platform.system(),
                run["cores"],
            )

            with tqdm.tqdm(unit="CDS", desc="HMMSCAN") as t:

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

        with tqdm.tqdm(unit="HSP", desc="HMMALIGN") as t:

            def align_callback(tasks_done: int):
                t.update(tasks_done)

            HMMer.align_simple(hsps_to_align, align_callback)

        alignment_count = 0
        for gbk in gbks:
            for gene in gbk.genes:
                for hsp in gene.hsps:
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
    # TODO: legacy bins

    # mix

    if not run["no_mix"] and not run["query_bgc_path"]:
        bs_mix.calculate_distances_mix(run, gbks)

        DB.commit()

    # legacy_classify

    if run["legacy_classify"] and not run["query_bgc_path"]:
        bs_legacy_classify.calculate_distances_legacy_classify(run, gbks)

        DB.commit()

    # classify

    if run["classify"] and not run["query_bgc_path"]:
        classify_mode = run["classify"]
        bs_classify.calculate_distances_classify(run, gbks, classify_mode)

        DB.commit()

    # query

    if run["query_bgc_path"]:
        bs_query.calculate_distances_query(run, gbks)

        DB.commit()

    # FAMILIES
    # TODO: per cutoff

    logging.info("Generating families")

    for cutoff in run["gcf_cutoffs"]:
        logging.info(" -- Cutoff %s", cutoff)
        for connected_component in bs_network.get_connected_components(cutoff):
            logging.debug(
                "Found connected component with %d edges", len(connected_component)
            )

            regions_families = bs_families.generate_families(
                connected_component, cutoff
            )

            # save families to database
            bs_families.save_to_db(regions_families)

    DB.commit()

    DB.save_to_disk(run["db_path"])

    # OUTPUT GENERATION

    # if this wasn't set in scan or align, set it now
    if pfam_info is None:
        HMMer.init(run["pfam_path"], False)
        pfam_info = HMMer.get_pfam_info()
        HMMer.unload()

    # prepare output files
    legacy_prepare_output(run["output_dir"], pfam_info)

    # prepare output files per cutoff
    for cutoff in run["gcf_cutoffs"]:
        legacy_prepare_cutoff_output(run["output_dir"], run["label"], cutoff, gbks)

    mix_bin = bs_comparison.RecordPairGenerator("mix")

    mix_bgc_regions: list[bs_gbk.BGCRecord] = []

    for gbk in gbks:
        if gbk.region is not None:
            gbk_records = bs_gbk.bgc_record.get_sub_records(
                gbk.region, run["record_type"]
            )
            mix_bgc_regions.extend(gbk_records)

    mix_bin.add_records(mix_bgc_regions)

    # mix

    if not run["no_mix"] and not run["query_bgc_path"]:
        mix_bin = bs_comparison.RecordPairGenerator("Mix")
        mix_bin.add_records([gbk.region for gbk in gbks if gbk.region is not None])

        for cutoff in run["gcf_cutoffs"]:
            legacy_prepare_bin_output(run["output_dir"], run["label"], cutoff, mix_bin)
            legacy_generate_bin_output(run["output_dir"], run["label"], cutoff, mix_bin)

    # legacy_classify

    if run["legacy_classify"]:
        legacy_class_bins = bs_comparison.legacy_bin_generator(gbks)

        for bin in legacy_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                legacy_prepare_bin_output(run["output_dir"], run["label"], cutoff, bin)
                legacy_generate_bin_output(run["output_dir"], run["label"], cutoff, bin)

    # classify

    if run["classify"]:
        classify_mode = run["classify"]

        if run["legacy_weights"]:
            weights = "legacy_weights"
        else:
            weights = "mix"

        as_class_bins = bs_comparison.as_class_bin_generator(
            gbks, weights, classify_mode
        )

        for bin in as_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                legacy_prepare_bin_output(run["output_dir"], run["label"], cutoff, bin)
                legacy_generate_bin_output(run["output_dir"], run["label"], cutoff, bin)

    # query
    # TODO: implement once classes are implemented for query mode

    # if run["query_bgc_path"]:
    #     query_mix_bin = bs_comparison.RecordPairGenerator("Query")
    #     query_mix_bin.add_records([gbk.region for gbk in gbks if gbk.region is not None])

    #     for cutoff in run["gcf_cutoffs"]:
    #         legacy_prepare_bin_output(run["output_dir"], run["label"], cutoff, query_mix_bin)
    #         legacy_generate_bin_output(run["output_dir"], run["label"], cutoff, query_mix_bin)

    if run["profiling"]:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
