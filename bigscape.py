# from python
import sys
import os
import logging
from datetime import datetime
import platform
from pathlib import Path

# from other modules
from src.data import DB
from src.genbank import SOURCE_TYPE, BGCRecord, CDS, GBK
from src.hmm import HMMer, legacy_filter_overlap, HSP
from src.parameters import RunParameters, parse_cmd
from src.comparison import (
    generate_mix,
    legacy_bin_generator,
    create_bin_network_edges_alt as create_bin_network_edges,
)
from src.file_input import load_dataset_folder, get_mibig, get_pfam
from src.diagnostics import Profiler
from src.network import BSNetwork
from src.output import (
    legacy_prepare_output,
    legacy_prepare_cutoff_output,
    legacy_prepare_bin_output,
    legacy_generate_bin_output,
)


def load_data(run: RunParameters):
    """Performs data loading tasks for BiG-SCAPE. This will reconstruct a run up to and
    including HMMAlign if a database file was detected in the run directory

    Args:
        run (RunParameters): Run parameters object for this run
    """


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
        
    # get reference if either MIBiG version or user-made reference dir passed
    if run.input.mibig_version:
        mibig_version_dir = get_mibig(run.input.mibig_version, bigscape_dir)
        # mibig_gbks = load_dataset_folder(mibig_version_dir, SOURCE_TYPE.MIBIG)

    if run.input.reference_dir:
        reference_gbks = load_dataset_folder(
            run.input.reference_dir, SOURCE_TYPE.REFERENCE
        )

    gbks = load_dataset_folder(
        run.input.input_dir,
        SOURCE_TYPE.QUERY,
        run.input.input_mode,
        run.input.include_gbk,
        run.input.exclude_gbk,
        run.input.cds_overlap_cutoff,
    )

    # get pfam
    # TODO: update to fetch from argument only
    get_pfam(run.input.pfam_version, run.input.pfam_path)

    # Load datasets
    gbks = load_dataset_folder(
        run.input.input_dir,
        SOURCE_TYPE.QUERY,
        run.input.input_mode,
        run.input.include_gbk,
        run.input.exclude_gbk,
        run.input.cds_overlap_cutoff,
        run.legacy,
    )

    exec_time = datetime.now() - start_time
    logging.info(
        "loaded %d gbks at %f seconds", len(gbks), exec_time.total_seconds()
    )

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

    if run.binning.mix:
        logging.info("Generating mix bin")

        mix_network = BSNetwork()
        all_regions: list[BGCRecord] = []
        for gbk in gbks:
            if gbk.region is not None:
                all_regions.append(gbk.region)
                mix_network.add_node(gbk.region)

        mix_bin = generate_mix(all_regions)

        legacy_prepare_bin_output(run.output.output_dir, run.label, 0.3, mix_bin)

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

        mix_network.generate_families_cutoff("dist", 0.3)

        # Output

        legacy_generate_bin_output(
            run.output.output_dir, run.label, 0.3, mix_bin, mix_network
        )

        mix_network.write_graphml(run.output.output_dir / Path("network_mix.graphml"))
        mix_network.write_edgelist_tsv(run.output.output_dir / Path("network_mix.tsv"))

        mix_network.export_distances_to_db()

    # networking - bins

    if not run.binning.legacy_no_classify:
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
            if run.binning.mix:
                continue

            bin_network.export_distances_to_db()

    # last dump of the database
    DB.save_to_disk(run.output.db_path)

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
