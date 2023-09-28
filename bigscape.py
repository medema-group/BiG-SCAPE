# from python
import sys
import os
import logging
from datetime import datetime
import platform
from pathlib import Path

# from other modules
from src.data import DB
from src.hmm import HMMer, legacy_filter_overlap
from src.parameters import RunParameters, parse_cmd
from src.diagnostics import Profiler
from src.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)


import src.file_input as bs_files
import src.data as bs_data
import src.enums as bs_enums
import src.comparison as bs_comparison
import src.network.network as bs_network
import src.network.families as bs_families
import src.distances.mix as bs_mix
import src.distances.query as bs_query


if __name__ == "__main__":
    bigscape_dir = Path(os.path.dirname(os.path.abspath(__file__)))

    # parsing needs to come first because we need it in setting up the logging
    run: RunParameters = parse_cmd(sys.argv[1:])

    start_time = run.start()

    # only now we can use logging.info etc to log stuff otherwise things get weird
    # initializing the logger and logger file also happens here
    run.validate()

    if run.legacy:
        logging.info("Using legacy mode")

    if not HMMer.are_profiles_pressed(run.input.pfam_path):
        logging.warning("HMM files were not pressed!")
        HMMer.press(run.input.pfam_path)

    # start profiler
    if run.diagnostics.profiling:
        profiler = Profiler(run.output.profile_path)
        profiler.start()

    # INPUT - load data

    gbks = bs_files.load_gbks(run, bigscape_dir)

    # get fist task
    run_state = bs_data.find_minimum_task(gbks)

    # nothing to do!
    if run_state == bs_enums.TASK.NOTHING_TO_DO:
        logging.info("Nothing to do!")
        exit(0)

    logging.info("First task: %s", run_state)

    # we will need this later. this can be set in hmmscan, hmmalign or not at all
    pfam_info = None

    # HMMER - Search/scan
    if run_state == bs_enums.TASK.HMM_SCAN:
        HMMer.init(run.input.pfam_path)

        # first opportunity to set this is here
        pfam_info = HMMer.get_pfam_info()

        # we certainly need to do some sort of scan, since the run_state is on hmm_scan. but we do
        # not need to do everything. Find the CDS that need scanning
        cds_to_scan = bs_data.get_cds_to_scan(gbks)

        logging.info("Scanning %d CDS", len(cds_to_scan))

        def callback(tasks_done):
            percentage = int(tasks_done / len(cds_to_scan) * 100)
            logging.info("%d/%d (%d%%)", tasks_done, len(cds_to_scan), percentage)

        if platform.system() == "Darwin":
            logging.warning("Running on mac-OS: hmmsearch_simple single threaded")
            HMMer.hmmsearch_simple(cds_to_scan, 1)
        else:
            logging.debug(
                "Running on %s: hmmsearch_multiprocess with %d cores",
                platform.system(),
                run.cores,
            )

        # TODO: the overlap filtering in this function does not seem to work
        HMMer.hmmsearch_multiprocess(
            cds_to_scan,
            domain_overlap_cutoff=run.hmmer.domain_overlap_cutoff,
            cores=run.cores,
            callback=callback,
        )

        # TODO: move, or remove after the add_hsp_overlap function is fixed (if it is broken
        # in the first place)
        # this sorts all CDS and then filters them using the old filtering system, which
        # is less efficient than the flitering using the CDS.add_hsp_overlap_filter
        # method. however, that method seems to be broken somehow
        all_hsps = []
        for gbk in gbks:
            for cds in gbk.genes:
                cds.hsps = sorted(cds.hsps)
                all_hsps.extend(
                    [hsp.domain for hsp in legacy_filter_overlap(cds.hsps, 0.1)]
                )

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

        all_hsps = []
        for gbk in gbks:
            for cds in gbk.genes:
                all_hsps.extend(cds.hsps)

        logging.info("%d hsps found in this run", len(all_hsps))

    # HMMER - Align
    if run_state == bs_enums.TASK.HMM_ALIGN:
        HMMer.init(run.input.pfam_path, False)

        # if this wasn't set before, set it now
        if pfam_info is None:
            pfam_info = HMMer.get_pfam_info()

        HMMer.align_simple(bs_data.get_hsp_to_align(gbks))

        alignment_count = 0
        for gbk in gbks:
            for gene in gbk.genes:
                for hsp in gene.hsps:
                    if hsp.alignment is None:
                        continue
                    hsp.alignment.save(False)

        logging.info("%d alignments", alignment_count)

        exec_time = datetime.now() - start_time
        logging.info("align done at %f seconds", exec_time.total_seconds())

        DB.commit()
        DB.save_to_disk(run.output.db_path)

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

    if not run.binning.no_mix and not run.binning.query_bgc_path:
        bs_mix.calculate_distances_mix(run, gbks)

        DB.commit()

    # query

    if run.binning.query_bgc_path:
        bs_query.calculate_distances_query(run, gbks)

        DB.commit()

    # FAMILIES
    # TODO: per cutoff

    logging.info("Generating families")

    for connected_component in bs_network.get_connected_components(0.3):
        logging.debug(
            "Found connected component with %d edges", len(connected_component)
        )

        regions_families = bs_families.generate_families(connected_component)

        # save families to database
        bs_families.save_to_db(regions_families)

    DB.commit()

    DB.save_to_disk(run.output.db_path)

    # OUTPUT

    # if this wasn't set in scan or align, set it now
    if pfam_info is None:
        HMMer.init(run.input.pfam_path, False)
        pfam_info = HMMer.get_pfam_info()
        HMMer.unload()

    # prepare output files
    legacy_prepare_output(run.output.output_dir, pfam_info)

    # work per cutoff
    legacy_prepare_cutoff_output(run.output.output_dir, run.label, 0.3, gbks)

    # TODO: I don't think the bins make much sense anymore
    # see if we can refactor this
    # TODO: per cutoff
    mix_bin = bs_comparison.RecordPairGenerator("mix")
    mix_bin.add_records([gbk.region for gbk in gbks if gbk.region is not None])

    legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

    legacy_generate_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
