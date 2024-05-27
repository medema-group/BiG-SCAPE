"""Contains functions for running pyhmmer's hmmscan workflow."""

# from python
from datetime import datetime
import logging
import platform
import tqdm
from typing import Any

# from other modules
from big_scape.data import DB
from big_scape.hmm import HMMer, HSP
import big_scape.data as bs_data


def run_hmmscan(run: dict[str, Any], gbks: list[Any], start_time: Any) -> None:
    """Runs the hmmscan workflow

    Args:
        run (dict): run parameters
        gbks (list[GBK]): list of GBK objects
        start_time (datetime): start time of the run
    """

    # HMMer.init(run["pfam_path"])

    # # first opportunity to set this is here
    # pfam_info = HMMer.get_pfam_info()

    # we certainly need to do some sort of scan, since the run_state is on hmm_scan. but we do
    # not need to do everything. Find the CDS that need scanning
    cds_to_scan = bs_data.get_cds_to_scan(gbks)

    logging.info("Scanning %d CDS", len(cds_to_scan))

    if platform.system() == "Darwin":
        logging.debug(
            "Running on %s: hmmsearch_simple with %d cores",
            platform.system(),
            run["cores"],
        )
        with tqdm.tqdm(unit="CDS", total=len(HMMer.profiles), desc="HMMSCAN") as t:
            HMMer.hmmsearch_simple(
                cds_to_scan,
                domain_overlap_cutoff=run["domain_overlap_cutoff"],
                cores=run["cores"],
                callback=lambda x, y: t.update(),
            )

    else:
        logging.debug(
            "Running on %s: hmmsearch_multiprocess with %d cores",
            platform.system(),
            run["cores"],
        )

        with tqdm.tqdm(unit="CDS", total=len(cds_to_scan), desc="HMMSCAN") as t:

            def callback(tasks_done):
                t.update(tasks_done)

            # TODO: the overlap filtering in this function does not seem to work
            HMMer.hmmsearch_multiprocess(
                cds_to_scan,
                domain_overlap_cutoff=run["domain_overlap_cutoff"],
                cores=run["cores"],
                callback=callback,
            )

    # sort hsps on starting position within the cds to ensure correct ordering during
    # lcs/extend and distance calculations
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

    all_hsps: list[HSP] = []
    for gbk in gbks:
        for cds in gbk.genes:
            all_hsps.extend(cds.hsps)

    logging.info("%d hsps found in this run", len(all_hsps))
