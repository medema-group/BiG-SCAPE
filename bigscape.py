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
from src.hmm import HMMer
from src.parameters import parse_cmd
from src.comparison import generate_mix
from src.comparison.legacy_extend import expand_glocal
from src.diagnostics import Profiler
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair
from src.network import BSNetwork

if __name__ == "__main__":
    # parsing needs to come first because we need it in setting up the logging
    run = parse_cmd(sys.argv[1:])

    # only now we can use logging.info etc to log stuff otherwise things get weird
    # initializing the logger and logger file also happens here
    run.validate()

    start_time = datetime.now()

    # start profiler
    if run.diagnostics.profiling:
        profiler = Profiler(run.output.profile_path)
        profiler.start()

    # start DB
    DB.create_in_mem()

    gbks = load_dataset_folder(
        run.input.input_dir,
        SOURCE_TYPE.QUERY,
        run.input.input_mode,
        run.input.include_gbk,
        run.input.exclude_gbk,
        run.input.cds_overlap_cutoff,
    )

    exec_time = datetime.now() - start_time
    logging.info("loaded %d gbks at %f seconds", len(gbks), exec_time.total_seconds())

    all_cds: list[CDS] = []
    for gbk in gbks:
        gbk.save_all()
        all_cds.extend(gbk.genes)

    logging.info("loaded %d cds total", len(all_cds))

    HMMer.init(run.input.pfam_path)

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
        HMMer.hmmsearch_multiprocess(all_cds, cores=run.cores, callback=callback)

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

    network = BSNetwork()
    all_regions: list[BGCRecord] = []
    for gbk in gbks:
        if gbk.region is not None:
            all_regions.append(gbk.region)
            network.add_node(gbk.region)

    mix_bin = generate_mix(all_regions)

    logging.info("Generated mix bin: %s", mix_bin)

    for pair in mix_bin.pairs():
        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            continue

        logging.debug("JC: %f", jaccard)

        pair.find_lcs()

        pair.comparable_region.log_comparable_region("LCS")

        expand_glocal(pair.comparable_region)

        pair.comparable_region.log_comparable_region("GLOCAL")

        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            continue

        adjacency = calc_ai_pair(pair)
        dss = calc_dss_pair(pair)

        # mix
        distance = 1 - (0.2 * jaccard) - (0.75 * adjacency) - (0.05 * dss)

        logging.debug(
            "JC: %f, AI: %f, DSS: %f, SCORE: %f", jaccard, adjacency, dss, distance
        )

        network.add_edge(pair, jc=jaccard, ai=adjacency, dss=dss, dist=distance)

    network.generate_families_cutoff("dist", 0.3)

    network.write_graphml(run.output.output_dir / Path("network.graphml"))

    if run.diagnostics.profiling:
        profiler.stop()

    exec_time = datetime.now() - start_time
    logging.info("All tasks done at %f seconds", exec_time.total_seconds())
