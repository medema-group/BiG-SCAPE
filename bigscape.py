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
from src.comparison import (
    RecordPairGenerator,
    RecordPairGeneratorQueryRef,
    PartialRecordPairGenerator,
    legacy_bin_generator,
    generate_edges,
    save_edge_to_db,
)
from src.file_input import load_dataset_folder, get_mibig
from src.diagnostics import Profiler
from src.network import BSNetwork, sim_matrix_from_edge_list, aff_sim_matrix
from src.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)

import src.data as bs_data
import src.enums as bs_enums


def load_gbks(run: RunParameters) -> list[GBK]:
    # counterintuitive, but we need to load the gbks first to see if there are any differences with
    # the data in the database

    print(sys.argv[0])

    input_gbks = []

    if run.binning.query_bgc_path:
        query_bgc_gbk = GBK(run.binning.query_bgc_path, bs_enums.SOURCE_TYPE.QUERY)
        input_gbks.append(query_bgc_gbk)

        gbks = load_dataset_folder(
            run.input.input_dir,
            bs_enums.SOURCE_TYPE.REFERENCE,
            run.input.input_mode,
            run.input.include_gbk,
            run.input.exclude_gbk,
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

        mix_bin = RecordPairGenerator("mix")
        mix_bin.add_records(mix_bgc_regions)

        logging.info(mix_bin)

        # check if there are existing distances
        missing_edge_bin = PartialRecordPairGenerator(mix_bin)
        num_missing_edges = missing_edge_bin.num_pairs()

        if num_missing_edges > 0:
            logging.info("Calculating distances for %d pairs", num_missing_edges)

            def callback(done_pairs):
                if mix_bin.num_pairs() > 10:
                    mod = round(mix_bin.num_pairs() / 10)
                else:
                    mod = 1

                if done_pairs % mod == 0:
                    logging.info(
                        f"{done_pairs}/{mix_bin.num_pairs()} ({done_pairs/mix_bin.num_pairs():.2%})"
                    )

            mix_edges = generate_edges(
                missing_edge_bin, run.comparison.alignment_mode, run.cores, callback
            )

            num_edges = 0
            for edge in mix_edges:
                num_edges += 1
                save_edge_to_db(edge)
                if num_edges == 100:
                    break

            logging.info("Generated %d edges", num_edges)

            DB.commit()
            DB.save_to_disk(run.output.db_path)

    exit(0)

    # NETWORKING

    # mix

    if not run.binning.no_mix and not run.binning.query_bgc_path:
        logging.info("Generating mix connected component networks")

        # labels will start at 0 and go up to the number of families. we have to
        # keep track of the last number we used, so we don't overlap numbers

        last_center = 0

        for connected_component in BSNetwork.generate_connected_component_networks(0.3):
            logging.info(
                "Found connected component with %d edges", len(connected_component)
            )

            logging.debug("Performing AP on connected component")

            distance_matrix = sim_matrix_from_edge_list(connected_component)

            labels, centers = aff_sim_matrix(distance_matrix)

            for label in labels:
                if label == -1:
                    continue

                family = label + last_center

            last_center = max(label)

            logging.debug("Found %d clusters", len(centers))

            # save families to database

        # mix_network.generate_families_cutoff("dist", 0.3)

        # # TODO: cull ref singletons again, for each family subgraph

        # # Output
        # legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

        # legacy_generate_bin_output(
        #     run.output.output_dir, run.label, 0.3, mix_bin, mix_network
        # )

        # mix_network.write_graphml(run.output.output_dir / Path("network_mix.graphml"))
        # mix_network.write_edgelist_tsv(run.output.output_dir / Path("network_mix.tsv"))

        # mix_network.export_distances_to_db()

    exit(0)

    # networking - query

    if run.binning.query_bgc_path:
        logging.info("Generating query BGC mode bin")

        query_bgc_network = BSNetwork()
        query_bgcs_records: list[BGCRecord] = []

        for gbk in gbks:
            if gbk.region is not None:
                query_bgcs_records.append(gbk.region)
                query_bgc_network.add_node(gbk.region)

        # all query <-> query and query <-> ref pairs
        query_bgc_pairs_queryref = RecordPairGeneratorQueryRef("Query_Ref")
        query_bgc_pairs_queryref.add_records(query_bgcs_records)

        #   create edges (query <-> query & ref edges)
        create_bin_network_edges(
            query_bgc_pairs_queryref,
            query_bgc_network,
            run.comparison.alignment_mode,
            run.cores,
            callback,
        )

        # create firstly all con_ref <-> sin_ref pairs and add edges in network
        query_bgc_pairs_conrefsinref = RecordPairGeneratorConRefSinRef(
            "ConRef_SinRef", query_bgc_network
        )
        query_bgc_pairs_conrefsinref.add_records(query_bgcs_records)
        create_bin_network_edges(
            query_bgc_pairs_conrefsinref,
            query_bgc_network,
            run.comparison.alignment_mode,
            run.cores,
            callback,
        )

        # get newly connected ref nodes, and recalculate nodes as to only compare new connected
        # refs to new singletons, to avoid double calculations
        query_bgc_pairs_conrefsinref.recalculate_nodes()

        old_nr_of_edges = 0
        new_nr_of_edges = query_bgc_network.graph.number_of_edges()

        # propagate network with new connected ref <-> new singletons until no new edges are added
        while new_nr_of_edges > old_nr_of_edges:
            old_nr_of_edges = new_nr_of_edges

            create_bin_network_edges(
                query_bgc_pairs_conrefsinref,
                query_bgc_network,
                run.comparison.alignment_mode,
                run.cores,
                callback,
            )

            query_bgc_pairs_conrefsinref.recalculate_nodes()
            new_nr_of_edges = query_bgc_network.graph.number_of_edges()

        #   cull leftover ref singletons
        query_bgc_network.cull_singletons(
            node_types=[bs_enums.SOURCE_TYPE.REFERENCE, bs_enums.SOURCE_TYPE.MIBIG]
        )

        all_current_nodes = query_bgc_network.get_nodes()

        # make pairs and add missing edges with all-vs-all pain_generator
        query_bgc_pairs_conrefconref = RecordPairGenerator("ConRef_ConRef")
        query_bgc_pairs_conrefsinref.add_records(all_current_nodes)
        create_bin_network_edges(
            query_bgc_pairs_conrefconref,
            query_bgc_network,
            run.comparison.alignment_mode,
            run.cores,
            callback,
        )

        # generate families
        # TODO: cull ref singletons again, for each family subgraph
        query_bgc_network.generate_families_cutoff("dist", 0.3)

        final_query_bgc_nodes = query_bgc_network.get_nodes()
        mock_query_bgc_pairs_conrefconref = RecordPairGenerator("Query_BGC")
        query_bgc_pairs_conrefsinref.add_records(final_query_bgc_nodes)

        # Output
        # TODO: check whether query_bgc_pairs_conrefsinref bin works here
        legacy_prepare_bin_output(
            run.output.output_dir, run.label, 0.3, query_bgc_pairs_conrefsinref
        )

        legacy_generate_bin_output(
            run.output.output_dir,
            run.label,
            0.3,
            query_bgc_pairs_conrefsinref,
            query_bgc_network,
        )

        query_bgc_network.write_graphml(
            run.output.output_dir / Path("network_mix.graphml")
        )
        query_bgc_network.write_edgelist_tsv(
            run.output.output_dir / Path("network_mix.tsv")
        )

        query_bgc_network.export_distances_to_db()

    # networking - bins

    if run.binning.legacy_classify:
        logging.info("Generating legacy bins")
        for bin in legacy_bin_generator(gbks):
            if bin.num_pairs() == 0:
                logging.info("Bin %s has no pairs. Skipping...", bin.label)

            bin_network = BSNetwork()
            for record in bin.source_records:
                bin_network.add_node(record)

            legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, bin)

            logging.info(bin)

            def callback(done_pairs):
                if bin.num_pairs() > 10:
                    mod = round(bin.num_pairs() / 10)
                else:
                    mod = 1

                if done_pairs % mod == 0:
                    logging.info(
                        f"{done_pairs}/{bin.num_pairs()} ({done_pairs/bin.num_pairs():.2%})"
                    )

            create_bin_network_edges(
                bin, bin_network, run.comparison.alignment_mode, run.cores, callback
            )

            bin_network.generate_families_cutoff("dist", 0.3)

            # Output

            legacy_generate_bin_output(
                run.output.output_dir, run.label, 0.3, bin, bin_network
            )

            bin_graphml_filename = f"network_{bin.label}.graphml"
            bin_network.write_graphml(
                run.output.output_dir / Path(bin_graphml_filename)
            )

            bin_tsv_filename = f"network_{bin.label}.tsv"
            bin_network.write_edgelist_tsv(
                run.output.output_dir / Path(bin_tsv_filename)
            )

            # if running mix, distances are already in the database
            if run.binning.no_mix:
                continue

            bin_network.export_distances_to_db()

    # last dump of the database
    DB.save_to_disk(run.output.db_path)

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
