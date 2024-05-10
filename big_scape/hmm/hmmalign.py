"""Contains functions for running pyhmmer's hmmalign workflow."""

# from python
from typing import Any
from datetime import datetime
import logging
import tqdm

# from other modules
from big_scape.data import DB
from big_scape.hmm import HMMer
import big_scape.data as bs_data


def run_hmmalign(run: dict[str, Any], gbks: list[Any], start_time: Any) -> None:
    """Runs the hmmalign workflow

    Args:
        run (dict): run parameters
        gbks (list[GBK]): list of GBK objects
        start_time (datetime): start time of the run
    """

    hsps_to_align = bs_data.get_hsp_to_align(gbks)

    with tqdm.tqdm(unit="HSP", total=len(hsps_to_align), desc="HMMALIGN") as t:

        def align_callback(tasks_done: int):
            t.update(tasks_done)

        HMMer.align_simple(hsps_to_align, align_callback)

    alignment_count = 0

    for hsp in hsps_to_align:
        if hsp.alignment is None:
            continue
        hsp.alignment.save(False)
        alignment_count += 1

    logging.info("%d alignments", alignment_count)

    exec_time = datetime.now() - start_time
    logging.info("align done at %f seconds", exec_time.total_seconds())

    DB.commit()
    DB.save_to_disk(run["db_path"])

    exec_time = datetime.now() - start_time
    logging.info("DB: HSP alignment save done at %f seconds", exec_time.total_seconds())

    HMMer.unload()
