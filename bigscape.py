# from python
import sys
import os
import logging
from datetime import datetime
import platform
from pathlib import Path

# from other modules
from src.data import DB
from src.genbank import BGCRecord, CDS, GBK
from src.hmm import HMMer, legacy_filter_overlap, HSP
from src.parameters import RunParameters, parse_cmd
from src.file_input import load_dataset_folder, get_mibig
from src.diagnostics import Profiler
from src.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)


import src.data as bs_data
import src.enums as bs_enums
import src.comparison as bs_comparison
import src.network.network as bs_network
import src.network.families as bs_families


def load_gbks(run: RunParameters) -> list[GBK]:
    # counterintuitive, but we need to load the gbks first to see if there are any differences with
    # the data in the database

    print(sys.argv[0])

    input_gbks = []

    if run.binning.query_bgc_path:
        query_bgc_gbk = GBK(run.binning.query_bgc_path, bs_enums.SOURCE_TYPE.QUERY)
        input_gbks.append(query_bgc_gbk)

        query_bgc_stem = run.binning.query_bgc_path.stem

        # add the query bgc to the exclude list
        exclude_gbk = run.input.exclude_gbk + [query_bgc_stem]

        gbks = load_dataset_folder(
            run.input.input_dir,
            bs_enums.SOURCE_TYPE.REFERENCE,
            run.input.input_mode,
            run.input.include_gbk,
            exclude_gbk,
            run.input.cds_overlap_cutoff,
        )
        input_gbks.extend(gbks)

    else:
        gbks = load_dataset_folder(
            run.input.input_dir,
            bs_enums.SOURCE_TYPE.QUERY,
            run.input.input_mode,
            run.input.include_gbk,
            run.input.exclude_gbk,
            run.input.cds_overlap_cutoff,
        )
        input_gbks.extend(gbks)

    # get reference if either MIBiG version or user-made reference dir passed
    if run.input.mibig_version:
        mibig_version_dir = get_mibig(run.input.mibig_version, bigscape_dir)
        mibig_gbks = load_dataset_folder(mibig_version_dir, bs_enums.SOURCE_TYPE.MIBIG)
        input_gbks.extend(mibig_gbks)

    if run.input.reference_dir:
        reference_gbks = load_dataset_folder(
            run.input.reference_dir, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(reference_gbks)

    # find the minimum task set for these gbks
    # if there is no database, create a new one and load in all the input stuff
    if not run.output.db_path.exists():
        DB.create_in_mem()

        all_cds: list[CDS] = []
        for gbk in input_gbks:
            gbk.save_all()
            all_cds.extend(gbk.genes)

        return input_gbks

    DB.load_from_disk(run.output.db_path)
    task_state = bs_data.find_minimum_task(input_gbks)

    # if we are are not on the load_gbks task, we have all the data we need
    if task_state != bs_enums.TASK.LOAD_GBKS:
        logging.info("Loading existing run from disk...")

        gbks = GBK.load_all()

        for gbk in gbks:
            HSP.load_all(gbk.genes)

        return gbks

    # if we end up here, we are in some halfway state and need to load in the new data
    logging.info("Loading existing run from disk and adding new data...")
    missing_gbks = bs_data.get_missing_gbks(input_gbks)
    logging.info("Found %d missing gbks", len(missing_gbks))

    for gbk in missing_gbks:
        gbk.save_all()

    # still return the full set
    return input_gbks


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

    gbks = load_gbks(run)

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

    # PREPARE OUTPUT

    # if this wasn't set in scan or align, set it now
    if pfam_info is None:
        HMMer.init(run.input.pfam_path, False)
        pfam_info = HMMer.get_pfam_info()
        HMMer.unload()

    # prepare output files
    legacy_prepare_output(run.output.output_dir, pfam_info)

    # work per cutoff
    # TODO: currently only 0.3
    legacy_prepare_cutoff_output(run.output.output_dir, run.label, 0.3, gbks)

    # DISTANCE GENERATION

    # mix

    if not run.binning.no_mix and not run.binning.query_bgc_path:
        logging.info("Generating mix bin")

        mix_bgc_regions: list[BGCRecord] = []

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

            DB.commit()

    # query

    if run.binning.query_bgc_path:
        logging.info("Generating query BGC mode bin")

        query_bgcs_records: list[BGCRecord] = []

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
                run.comparison.alignment_mode,
                run.cores,
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
                        f"{done_pairs}/{num_pairs} ({done_pairs/num_pairs:.2%})"
                    )

            ref_edges = bs_comparison.generate_edges(
                ref_to_ref_bin, run.comparison.alignment_mode, run.cores, callback
            )

            num_edges = 0
            for edge in ref_edges:
                num_edges += 1
                bs_comparison.save_edge_to_db(edge, True)

            if num_edges == 0:
                break

            logging.info("Generated %d edges", num_edges)

        DB.commit()

    # NETWORKING

    # mix

    logging.info("Generating mix connected component networks")

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

    # mix

    if not run.binning.no_mix and not run.binning.query_bgc_path:
        legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

        legacy_generate_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

    if run.binning.query_bgc_path:
        legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, ref_to_ref_bin)

        legacy_generate_bin_output(
            run.output.output_dir, run.label, 0.3, ref_to_ref_bin
        )

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
