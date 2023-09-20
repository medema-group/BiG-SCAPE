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
    RecordPairGeneratorConRefSinRef,
    RecordPairGeneratorGivenNodeSinRef,
    generate_mix,
    legacy_bin_generator,
    create_bin_network_edges_alt as create_bin_network_edges,
    # Alt is an alternative multiprocessing implementation
)
from src.file_input import load_dataset_folder, get_mibig
from src.diagnostics import Profiler
from src.network import BSNetwork
from src.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)
from src.enums import SOURCE_TYPE

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

    if run.output.db_path.exists():
        logging.info("Loading existing run from disk...")

        DB.load_from_disk(run.output.db_path)

        gbks = GBK.load_all()

        for gbk in gbks:
            HSP.load_all(gbk.genes)

        HMMer.init(run.input.pfam_path)

        pfam_info = HMMer.get_pfam_info()

        HMMer.unload()

    else:
        # start DB
        DB.create_in_mem()

    if run.binning.query_bgc_path:
        query_bgc_gbk = GBK(run.binning.query_bgc_path, SOURCE_TYPE.QUERY)

        gbks = load_dataset_folder(
            run.input.input_dir,
            SOURCE_TYPE.REFERENCE,
            run.input.input_mode,
            run.input.include_gbk,
            run.input.exclude_gbk,
            run.input.cds_overlap_cutoff,
        )

        gbks.append(query_bgc_gbk)

    else:
        gbks = load_dataset_folder(
            run.input.input_dir,
            SOURCE_TYPE.QUERY,
            run.input.input_mode,
            run.input.include_gbk,
            run.input.exclude_gbk,
            run.input.cds_overlap_cutoff,
        )

    # get reference if either MIBiG version or user-made reference dir passed
    if run.input.mibig_version:
        mibig_version_dir = get_mibig(run.input.mibig_version, bigscape_dir)
        mibig_gbks = load_dataset_folder(mibig_version_dir, SOURCE_TYPE.MIBIG)
        gbks.extend(mibig_gbks)

    if run.input.reference_dir:
        reference_gbks = load_dataset_folder(
            run.input.reference_dir, SOURCE_TYPE.REFERENCE
        )
        gbks.extend(reference_gbks)

    all_cds: list[CDS] = []
    for gbk in gbks:
        gbk.save_all()
        all_cds.extend(gbk.genes)

    logging.info("loaded %d cds total", len(all_cds))

    # HMMER - Search

    HMMer.init(run.input.pfam_path)

    # we will need this information for the output later
    pfam_info = HMMer.get_pfam_info()

    def callback(tasks_done):
        percentage = int(tasks_done / len(all_cds) * 100)
        logging.info("%d/%d (%d%%)", tasks_done, len(all_cds), percentage)

    if platform.system() == "Darwin":
        logging.warning("Running on mac-OS: hmmsearch_simple single threaded")
        HMMer.hmmsearch_simple(all_cds, 1)
    else:
        logging.debug(
            "Running on %s: hmmsearch_multiprocess with %d cores",
            platform.system(),
            run.cores,
        )

        # TODO: the overlap filtering in this function does not seem to work
        HMMer.hmmsearch_multiprocess(
            all_cds,
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

        all_hsps = []
        for cds in all_cds:
            all_hsps.extend(cds.hsps)

        logging.info("%d hsps", len(all_hsps))

        exec_time = datetime.now() - start_time
        logging.info("scan done at %f seconds", exec_time.total_seconds())

        HMMer.unload()

        # save hsps to database
        for new_hsp in all_hsps:
            new_hsp.save(False)
        DB.commit()

        exec_time = datetime.now() - start_time
        logging.info("DB: HSP save done at %f seconds", exec_time.total_seconds())

        # HMMER - Align

        HMMer.init(run.input.pfam_path, False)

        HMMer.align_simple(all_hsps)

        all_alignments = list()
        for cds in all_cds:
            for hsp in cds.hsps:
                if hsp.alignment is None:
                    continue
                all_alignments.append(hsp.alignment)

        logging.info("%d alignments", len(all_alignments))

        exec_time = datetime.now() - start_time
        logging.info("align done at %f seconds", exec_time.total_seconds())

        HMMer.unload()

        for hsp_alignment in all_alignments:
            hsp_alignment.save(False)
        DB.commit()

        exec_time = datetime.now() - start_time
        logging.info(
            "DB: HSP alignment save done at %f seconds", exec_time.total_seconds()
        )

        DB.save_to_disk(run.output.db_path)

    # prepare output files
    legacy_prepare_output(run.output.output_dir, pfam_info)

    # work per cutoff
    # TODO: currently only 0.3
    legacy_prepare_cutoff_output(run.output.output_dir, run.label, 0.3, gbks)

    # networking - mix

    if not run.binning.no_mix:
        logging.info("Generating mix bin")

        mix_bgc_records: list[BGCRecord] = []
        mix_network = BSNetwork()

        for gbk in gbks:
            if gbk.region is not None:
                mix_bgc_records.append(gbk.region)
                mix_network.add_node(gbk.region)

        # needs to be generalized
        mix_bin = generate_mix(mix_bgc_records)

        logging.info(mix_bin)

        def callback(done_pairs):
            if mix_bin.num_pairs() > 10:
                mod = round(mix_bin.num_pairs() / 10)
            else:
                mod = 1

            if done_pairs % mod == 0:
                logging.info(
                    f"{done_pairs}/{mix_bin.num_pairs()} ({done_pairs/mix_bin.num_pairs():.2%})"
                )

        create_bin_network_edges(
            mix_bin, mix_network, run.comparison.alignment_mode, run.cores, callback
        )

        mix_network.cull_singletons(
            node_types=[SOURCE_TYPE.MIBIG, SOURCE_TYPE.REFERENCE]
        )

        mix_network.generate_families_cutoff("dist", 0.3)

        # TODO: cull ref singletons again, for each family subgraph

        # Output
        legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

        legacy_generate_bin_output(
            run.output.output_dir, run.label, 0.3, mix_bin, mix_network
        )

        mix_network.write_graphml(run.output.output_dir / Path("network_mix.graphml"))
        mix_network.write_edgelist_tsv(run.output.output_dir / Path("network_mix.tsv"))

        mix_network.export_distances_to_db()

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
        query_bgc_pairs_queryref.add_bgcs(query_bgcs_records)

        #   create edges (query <-> query & ref edges)
        create_bin_network_edges(
            query_bgc_pairs_queryref,
            query_bgc_network,
            run.comparison.alignment_mode,
            run.cores,
            callback,
        )

        ref_singletons = query_bgc_network.get_singletons(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )

        ref_nodes = query_bgc_network.get_nodes(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )

        con_ref_nodes = [node for node in ref_nodes if node not in ref_singletons]

        # create all con_ref <-> sin_ref pairs and add edges in network
        query_bgc_pairs_conrefsinref = RecordPairGeneratorConRefSinRef(
            "ConRef_SinRef", query_bgc_network
        )
        query_bgc_pairs_conrefsinref.add_bgcs(ref_singletons)
        query_bgc_pairs_conrefsinref.add_bgcs(con_ref_nodes)
        create_bin_network_edges(
            query_bgc_pairs_conrefsinref,
            query_bgc_network,
            run.comparison.alignment_mode,
            run.cores,
            callback,
        )

        # get new_ref_connected, by checking which ref_singletons are no longer available
        new_ref_singletons = query_bgc_network.get_singletons(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )
        new_ref_connected = [
            node for node in ref_singletons if node not in new_ref_singletons
        ]

        old_nr_of_edges = 0
        new_nr_of_edges = query_bgc_network.graph.number_of_edges()

        while new_nr_of_edges > old_nr_of_edges:
            old_nr_of_edges = new_nr_of_edges
            ref_singletons = new_ref_singletons

            # create pairs and edges
            query_bgc_pairs_givensinref = RecordPairGeneratorGivenNodeSinRef(
                "GivenRef_SinRef", query_bgc_network, new_ref_connected
            )
            query_bgc_pairs_givensinref.add_bgcs(new_ref_singletons)
            query_bgc_pairs_givensinref.add_bgcs(new_ref_connected)
            create_bin_network_edges(
                query_bgc_pairs_givensinref,
                query_bgc_network,
                run.comparison.alignment_mode,
                run.cores,
                callback,
            )

            # stored last updated edges
            new_ref_singletons = query_bgc_network.get_singletons(
                node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
            )
            new_ref_connected = [
                node for node in ref_singletons if node not in new_ref_singletons
            ]

            # new nr of edges = network.graph.number_of_edges()
            new_nr_of_edges = query_bgc_network.graph.number_of_edges()

        #   cull leftover ref singletons
        query_bgc_network.cull_singletons(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )

        # make sure all con-ref to con-ref edges are added
        # get all ref nodes and sin-ref nodes, and con_ref nodes from there
        ref_nodes = query_bgc_network.get_nodes(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )
        singleton_ref_nodes = query_bgc_network.get_nodes(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )
        connected_ref_nodes = [
            node for node in ref_nodes if node not in singleton_ref_nodes
        ]

        # make pairs and add edges with all-vs-all pain_generator
        query_bgc_pairs_conrefconref = RecordPairGenerator("ConRef_ConRef")
        query_bgc_pairs_conrefsinref.add_bgcs(connected_ref_nodes)
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
        query_bgc_pairs_conrefsinref.add_bgcs(final_query_bgc_nodes)

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
