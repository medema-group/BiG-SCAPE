"""The highest-level module of BiG-SCAPE.
This module is responsible for running the BiG-SCAPE analysis.
"""

# from python
import sys
import os
import psutil
import signal
import logging
from datetime import datetime
from pathlib import Path

# from other modules
from big_scape.hmm import HMMer, run_hmmsearch, run_hmmalign
from big_scape.diagnostics import Profiler
from big_scape.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
    write_record_annotations_file,
    write_full_network_file,
)
from big_scape.utility import domain_includelist_filter, class_filter, category_filter


import big_scape.file_input as bs_files

# import big_scape.genbank as bs_gbk
from big_scape.utility.version import get_bigscape_version
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
    # this is, lightly and delicately put, *not a great way to get the version*
    # but after some research, it seems no one had the idea to make this easy

    bigscape_version = get_bigscape_version()

    logging.info(
        "Starting BiG-SCAPE %s %s run on %s level with %s alignment and %s weights",
        bigscape_version,
        run["mode"].value,
        run["record_type"].value,
        run["alignment_mode"].value,
        "legacy" if run["legacy_weights"] else "mix",
    )

    # root directory
    bigscape_dir = Path(os.path.dirname(os.path.abspath(__file__)))

    main_pid = os.getpid()

    # capture ctrl-c for database saving and other cleanup
    def signal_handler(sig, frame):
        nonlocal main_pid

        # return if not main process
        if os.getpid() != main_pid:
            return
        if bs_data.DB.opened():
            logging.warning("User requested SIGINT. Saving database and exiting")
            logging.warning("Press ctrl+c again to force exit")

            # kill all child processes
            for child in psutil.Process().children(recursive=True):
                child.kill()

            # reset to default handler in case we get another SIGINT
            signal.signal(signal.SIGINT, signal.SIG_DFL)

            bs_data.DB.commit()
            bs_data.DB.save_to_disk(run["db_path"])

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

    # INPUT - create an in memory bs_data.DB or load from disk
    if not run["db_path"].exists():
        if run["disk_only"]:
            logging.info("Creating on disk database")
            bs_data.DB.create_on_disk(run["db_path"])
        else:
            logging.info("Creating in memory database")
            bs_data.DB.create_in_mem()
    else:
        logging.info("Loading database from disk")
        bs_data.DB.load_from_disk(run["db_path"])
        bs_data.DB.check_config_hash()

    # initialize run in database
    bs_data.DB.init_run(run)

    # INPUT - load data
    gbks = bs_files.load_gbks(run, bigscape_dir)

    # INPUT - get all working BGC records
    if run["query_bgc_path"]:
        all_bgc_records, query_record = bs_files.get_all_bgc_records_query(run, gbks)
    else:
        all_bgc_records = bs_files.get_all_bgc_records(run, gbks)

    logging.info("Continuing with %s BGC records", len(all_bgc_records))

    # INCLUDE/EXCLUDE CATEGORIES
    if run["include_categories"] or run["exclude_categories"]:
        all_bgc_records = category_filter(run, all_bgc_records)

    # INCLUDE/EXCLUDE CLASSES
    if run["include_classes"] or run["exclude_classes"]:
        all_bgc_records = class_filter(run, all_bgc_records)

    # get fist task
    run_state = bs_data.find_minimum_task(gbks)

    logging.info("First task: %s", run_state)

    # we will need this later. this can be set in hmmsearch, hmmalign or not at all
    pfam_info = None

    # HMMER - Search/scan
    if run_state == bs_enums.TASK.HMM_SCAN:
        HMMer.init(run["pfam_path"])

        # first opportunity to set this is here
        pfam_info = HMMer.get_pfam_info()

        run_hmmsearch(run, gbks, start_time)

        # set new run state
        run_state = bs_data.find_minimum_task(gbks)
        logging.info("Next task: %s", run_state)

    # HMMER - Align
    if run_state == bs_enums.TASK.HMM_ALIGN:
        HMMer.init(run["pfam_path"], False)

        # if this wasn't set before, set it now
        if pfam_info is None:
            pfam_info = HMMer.get_pfam_info()

        run_hmmalign(run, gbks, start_time)

        # set new run state
        run_state = bs_data.find_minimum_task(gbks)
        logging.info("Next task: %s", run_state)

    # TODO: idea: use sqlite to set distances of 1.0 for all pairs that have no domains
    # in common

    # DOMAIN INCLUSION LIST FILTER
    if run["domain_includelist_all"] or run["domain_includelist_any"]:
        all_bgc_records = domain_includelist_filter(run, all_bgc_records)

        logging.info(
            "Continuing with %i BGC records after domain_includelist filtering",
            len(all_bgc_records),
        )

    # DISTANCE GENERATION

    # mix

    if run["mix"] and not run["query_bgc_path"]:
        bs_mix.calculate_distances_mix(run, all_bgc_records)

        bs_data.DB.commit()

    # legacy_classify

    if run["classify"] == bs_enums.CLASSIFY_MODE.LEGACY and not run["query_bgc_path"]:
        bs_legacy_classify.calculate_distances_legacy_classify(run, all_bgc_records)

        bs_data.DB.commit()

    # classify

    if (
        run["classify"]
        and run["classify"] != bs_enums.CLASSIFY_MODE.LEGACY
        and not run["query_bgc_path"]
    ):
        bs_classify.calculate_distances_classify(run, all_bgc_records)

        bs_data.DB.commit()

    # query

    if run["query_bgc_path"]:
        query_bin = bs_query.calculate_distances_query(
            run, all_bgc_records, query_record
        )

        bs_data.DB.commit()

    bs_data.DB.save_to_disk(run["db_path"])

    # FAMILY GENERATION
    logging.info("Generating families")

    # mix

    if run["mix"] and not run["query_bgc_path"]:
        bs_families.run_family_assignments(
            run, bs_comparison.generate_mix_bin, all_bgc_records
        )

    # legacy_classify

    if run["classify"] == bs_enums.CLASSIFY_MODE.LEGACY and not run["query_bgc_path"]:
        bs_families.run_family_assignments(
            run, bs_comparison.legacy_bin_generator, all_bgc_records
        )

    # classify

    if (
        run["classify"]
        and run["classify"] != bs_enums.CLASSIFY_MODE.LEGACY
        and not run["query_bgc_path"]
    ):
        bs_families.run_family_assignments(
            run, bs_comparison.as_class_bin_generator, all_bgc_records
        )

    # query BGC mode
    if run["query_bgc_path"]:
        cc_cutoff = bs_families.run_family_assignments_query(
            run, query_bin, query_record
        )

    bs_data.DB.save_to_disk(run["db_path"])

    # OUTPUT GENERATION

    # if this wasn't set in scan or align, set it now
    if pfam_info is None:
        HMMer.init(run["pfam_path"], False)
        pfam_info = HMMer.get_pfam_info()
        HMMer.unload()

    logging.info("Generating GCF trees and output files")

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

    if run["mix"] and not run["query_bgc_path"]:
        mix_bin = bs_comparison.generate_mix_bin(all_bgc_records, run)

        for cutoff in run["gcf_cutoffs"]:
            mix_bin.cull_singletons(cutoff, run["run_id"])
            if len(mix_bin.record_ids) == 0:
                logging.info(
                    f"Network '{mix_bin.label}' with cutoff {cutoff} is empty after culling singletons"
                )
                continue
            legacy_prepare_bin_output(run, cutoff, mix_bin)
            legacy_generate_bin_output(run, cutoff, mix_bin)

    # legacy_classify

    if run["classify"] == bs_enums.CLASSIFY_MODE.LEGACY and not run["query_bgc_path"]:
        legacy_class_bins = bs_comparison.legacy_bin_generator(all_bgc_records, run)

        for bin in legacy_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                bin.cull_singletons(cutoff, run["run_id"])
                if len(bin.record_ids) == 0:
                    logging.info(
                        f"Network '{bin.label}' with cutoff {cutoff} is empty after culling singletons"
                    )
                    continue
                legacy_prepare_bin_output(run, cutoff, bin)
                legacy_generate_bin_output(run, cutoff, bin)

    # classify

    if (
        run["classify"]
        and run["classify"] != bs_enums.CLASSIFY_MODE.LEGACY
        and not run["query_bgc_path"]
    ):
        as_class_bins = bs_comparison.as_class_bin_generator(all_bgc_records, run)

        for bin in as_class_bins:
            for cutoff in run["gcf_cutoffs"]:
                bin.cull_singletons(cutoff, run["run_id"])
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
            if cutoff not in cc_cutoff.keys():
                continue
            query_connected_component = cc_cutoff[cutoff]
            query_records = bs_network.get_nodes_from_cc(
                query_connected_component, query_bin.source_records
            )

            query_bin = bs_comparison.ConnectedComponentPairGenerator(
                cc_cutoff[cutoff], query_bin.label
            )
            query_bin.add_records(
                [record for record in query_records if record is not None]
            )
            legacy_prepare_bin_output(run, cutoff, query_bin)
            legacy_generate_bin_output(run, cutoff, query_bin)

    if run["profiling"]:
        profiler.stop()

    end = datetime.now()
    exec_time = end - start_time

    bs_data.DB.set_run_end(run["run_id"], start_time, end)
    bs_data.DB.save_to_disk(run["db_path"], True)

    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
