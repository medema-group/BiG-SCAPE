# from python
import sys
import logging
from datetime import datetime
import platform
from pathlib import Path

# from other modules
from src.data import DB
from src.file_input import load_dataset_folder
from src.genbank import SOURCE_TYPE, BGCRecord, CDS
from src.hmm import HMMer, legacy_filter_overlap
from src.parameters import RunParameters, parse_cmd
from src.comparison import generate_mix, legacy_bin_generator, create_bin_network_edges
from src.diagnostics import Profiler
from src.network import BSNetwork
from src.output import generate_legacy_output


if __name__ == "__main__":
    # initialize run

    # parsing needs to come first because we need it in setting up the logging
    run: RunParameters = parse_cmd(sys.argv[1:])

    # only now we can use logging.info etc to log stuff otherwise things get weird
    # initializing the logger and logger file also happens here
    run.validate()

    if run.legacy:
        logging.info("Using legacy mode")

    start_time = datetime.now()

    # start profiler
    if run.diagnostics.profiling:
        profiler = Profiler(run.output.profile_path)
        profiler.start()

    # start DB
    DB.create_in_mem()

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
    logging.info("loaded %d gbks at %f seconds", len(gbks), exec_time.total_seconds())

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
        all_hsps = list(HMMer.hmmsearch_simple(all_cds, 1))
    else:
        logging.debug(
            "Running on %s: hmmsearch_multiprocess with %d cores",
            platform.system(),
            run.cores,
        )

        # if legacy is true, set cutoff to 1.1 for the domain filtering so we can use
        # legacy filtering later

        if run.legacy:
            domain_overlap_cutoff = 1.1
        else:
            domain_overlap_cutoff = run.hmmer.domain_overlap_cutoff
        HMMer.hmmsearch_multiprocess(
            all_cds,
            domain_overlap_cutoff=domain_overlap_cutoff,
            cores=run.cores,
        )

    # TODO: move, or remove after the add_hsp_overlap function is fixed (if it is broken
    # in the first place)
    for cds in all_cds:
        if run.legacy:
            cds.hsps = legacy_filter_overlap(cds.hsps, 0.1)

        # TODO: remove when sortedlists are removed
        # this converts the sortedlist used internally to regular lists.
        # for some reason, doing this beforehand really messes things up for reasons I don't
        # understand.
        cds.lock()

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
            all_alignments.append(hsp.alignment)

    logging.info("%d alignments", len(all_alignments))

    exec_time = datetime.now() - start_time
    logging.info("align done at %f seconds", exec_time.total_seconds())

    HMMer.unload()

    for hsp_alignment in all_alignments:
        hsp_alignment.save(False)
    DB.commit()

    exec_time = datetime.now() - start_time
    logging.info("DB: HSP alignment save done at %f seconds", exec_time.total_seconds())

    DB.save_to_disk(run.output.db_path)

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

        logging.info(mix_bin)

        create_bin_network_edges(mix_bin, mix_network, run.comparison.alignment_mode)

        mix_network.generate_families_cutoff("dist", 0.3)

        # Output

        mix_network.write_graphml(run.output.output_dir / Path("network_mix.graphml"))
        mix_network.write_edgelist_tsv(run.output.output_dir / Path("network_mix.tsv"))

        generate_legacy_output(
            run.output.output_dir, "test", [0.3], ["mix"], mix_network, gbks, pfam_info
        )

    # networking - bins

    if not run.binning.legacy_no_classify:
        logging.info("Generating legacy bins")
        for bin in legacy_bin_generator(gbks):
            if bin.num_pairs() == 0:
                logging.info("Bin %s has no pairs. Skipping...", bin.label)

            bin_network = BSNetwork()
            for record in bin.source_records:
                bin_network.add_node(record)

            logging.info(bin)

            create_bin_network_edges(bin, bin_network, run.comparison.alignment_mode)

            bin_network.generate_families_cutoff("dist", 0.3)

            # Output

            bin_graphml_filename = f"network_{bin.label}.graphml"
            bin_network.write_graphml(
                run.output.output_dir / Path(bin_graphml_filename)
            )

            bin_tsv_filename = f"network_{bin.label}.tsv"
            bin_network.write_edgelist_tsv(
                run.output.output_dir / Path(bin_tsv_filename)
            )

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
