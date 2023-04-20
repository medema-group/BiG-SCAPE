# from python
import sys
import logging
from datetime import datetime

# from other modules
from src.data import DB
from src.file_input import load_dataset_folder
from src.genbank import SOURCE_TYPE, BGCRecord, CDS
from src.hmm import HMMer, HSP
from src.parameters import parse_cmd
from src.comparison import generate_mix

if __name__ == "__main__":
    # parsing needs to come first because we need it in setting up the logging
    run = parse_cmd(sys.argv[1:])

    # logger
    # this tells the logger what the messages should look like
    # asctime = YYYY-MM-DD HH:MM:SS,fff
    # levelname = DEBUG/INFO/WARN/ERROR
    # message = whatever we pass, eg logging.debug("message")
    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")

    # get the built in logger
    root_logger = logging.getLogger()

    if run.diagnostics.verbose:
        root_logger.level = logging.DEBUG
    else:
        root_logger.level = logging.INFO

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)

    # only now we can use logging.info etc to log stuff otherwise things get weird
    run.validate()

    start_time = datetime.now()

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

    HMMer.init(run.input.pfam_path)

    all_cds: list[CDS] = []
    for gbk in gbks:
        gbk.save_all()
        all_cds.extend(gbk.genes)

    def callback(tasks_done):
        percentage = int(tasks_done / len(all_cds) * 100)
        logging.info("%d/%d (%d%%)", tasks_done, len(all_cds), percentage)

    HMMer.hmmsearch_multiprocess(all_cds)

    all_hsps: list[HSP] = []
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

    all_alignments = list(HMMer.align_simple(all_hsps))

    logging.info("%d alignments", len(all_alignments))

    exec_time = datetime.now() - start_time
    logging.info("align done at %f seconds", exec_time.total_seconds())

    for hsp_alignment in all_alignments:
        hsp_alignment.save(False)
    DB.commit()

    exec_time = datetime.now() - start_time
    logging.info("DB: HSP alignment save done at %f seconds", exec_time.total_seconds())

    DB.save_to_disk(run.output.db_path)

    all_regions: list[BGCRecord] = []
    for gbk in gbks:
        if gbk.region is not None:
            all_regions.append(gbk.region)

    mix_bin = generate_mix(all_regions)

    logging.info("Generated mix bin: %s", mix_bin)
