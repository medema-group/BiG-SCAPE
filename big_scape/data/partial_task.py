"""Contains code that scans the database for existing data and determines at which stage
of executing BiG-SCAPE the application last exited

Also contains functions to determine subsets of tasks that need to be done
"""

# from python
from __future__ import annotations
from typing import TYPE_CHECKING

# from dependencies
from sqlalchemy import select

# from other modules
from big_scape.data import DB
import big_scape.enums as bs_enums

# from circular imports
if TYPE_CHECKING:
    from big_scape.genbank.gbk import GBK, CDS
    from big_scape.hmm import HSP


def find_minimum_task(gbks: list[GBK]) -> bs_enums.TASK:
    """Finds the earliest bs_enums.TASK to start at based on data in the database

    If new data was added, this will always be the load_gbks bs_enums.TASK,
    otherwise, it tries to find the latest bs_enums.TASK with unfinished business

    Args:
        gbks (list[GBK]): input_gbks

    Returns:
        bs_enums.TASK: minimum task
    """
    input_data_state = get_input_data_state(gbks)

    # new data or mixed data
    if (
        input_data_state.value == bs_enums.INPUT_TASK.NEW_DATA.value
        or input_data_state.value == bs_enums.INPUT_TASK.MIXED_DATA.value
        or input_data_state.value == bs_enums.INPUT_TASK.NO_DATA.value
    ):
        # gbks from input need to be loaded into the in-memory database
        return bs_enums.TASK.SAVE_GBKS

    hmm_data_state = get_hmm_data_state()

    if hmm_data_state.value < bs_enums.HMM_TASK.NEED_ALIGN.value:
        return bs_enums.TASK.HMM_SCAN

    if hmm_data_state.value < bs_enums.HMM_TASK.ALL_ALIGNED.value:
        return bs_enums.TASK.HMM_ALIGN

    # if scan and align are finished, we are ready for comparison
    return bs_enums.TASK.COMPARISON


def get_input_data_state(gbks: list[GBK]) -> bs_enums.INPUT_TASK:
    """Returns the status of input data (gbks and regions) in the in-memory database

    Args:
        gbks (list[GBK]): input gbks

    Returns:
        bs_enums.INPUT_TASK: status of input data
    """
    gbk_count = DB.get_table_row_count("gbk")

    if gbk_count == 0:
        return bs_enums.INPUT_TASK.NO_DATA

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    gbk_table = DB.metadata.tables["gbk"]

    # get set of gbks in database
    db_gbk_rows = DB.execute(select(gbk_table.c.hash))
    db_gbk_hashes: set[str] = {row[0] for row in db_gbk_rows.all()}
    input_gbk_hashes: set[str] = {str(gbk.hash) for gbk in gbks}

    if db_gbk_hashes == input_gbk_hashes:
        return bs_enums.INPUT_TASK.SAME_DATA

    union = db_gbk_hashes & input_gbk_hashes

    # all new data
    if len(union) == 0:
        return bs_enums.INPUT_TASK.NEW_DATA

    # only partial data which is already in database
    if len(union) == len(input_gbk_hashes):
        return bs_enums.INPUT_TASK.PARTIAL_DATA

    # otherwise there is some new data, some old data is missing
    return bs_enums.INPUT_TASK.MIXED_DATA


def get_missing_gbks(gbks: list[GBK]) -> list[GBK]:
    """Find which GBKs are missing from the database and return them

    Args:
        gbks (list[GBK]): List of GBKs to check

    Returns:
        list[GBK]: List of GBKs that are missing from the database
    """
    # dictionary of gbk path to gbk object
    gbk_dict = {str(gbk.hash): gbk for gbk in gbks}

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    gbk_table = DB.metadata.tables["gbk"]

    # get set of gbks in database
    db_gbk_rows = DB.execute(select(gbk_table.c.hash))
    db_gbk_hashes: set[str] = {row[0] for row in db_gbk_rows.all()}

    missing_gbks = []

    for gbk_hash in gbk_dict:
        if gbk_hash not in db_gbk_hashes:
            missing_gbks.append(gbk_dict[gbk_hash])

    return missing_gbks


def get_hmm_data_state() -> bs_enums.HMM_TASK:
    """Returns the state of data hmm processing in the in-memory database

    Returns:
        bs_enums.HMM_TASK: state of hmm data
    """
    hsp_count = DB.get_table_row_count("hsp")

    if hsp_count == 0:
        return bs_enums.HMM_TASK.NO_DATA

    cds_count = DB.get_table_row_count("cds")
    scanned_cds_count = DB.get_table_row_count("scanned_cds")

    if cds_count > scanned_cds_count:
        return bs_enums.HMM_TASK.NEED_SCAN

    align_count = DB.get_table_row_count("hsp_alignment")

    if align_count < hsp_count:
        return bs_enums.HMM_TASK.NEED_ALIGN

    return bs_enums.HMM_TASK.ALL_ALIGNED


def get_cds_to_scan(gbks: list[GBK]) -> list[CDS]:
    """Find which cds within a list of GBK objects need to be scanned using hmmscan

    Args:
        gbks (list[GBK]): List of GBKs

    Returns:
        list[CDS]: List of CDS which were not scanned according to the data in scanned_cds table
    """
    hmm_state = get_hmm_data_state()

    if hmm_state == bs_enums.HMM_TASK.NO_DATA:
        cds_to_scan = []
        for gbk in gbks:
            cds_to_scan.extend(gbk.genes)
        return cds_to_scan

    cds_to_scan = []

    # get a list of database cds_ids that are present in the cds_scanned table

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    scanned_cds_table = DB.metadata.tables["scanned_cds"]
    select_query = select(scanned_cds_table.c.cds_id)
    scanned_cds_ids = set(DB.execute(select_query))

    for gbk in gbks:
        for gene in gbk.genes:
            if (gene._db_id,) not in scanned_cds_ids:
                cds_to_scan.append(gene)

    return cds_to_scan


def get_hsp_to_align(gbks: list[GBK]) -> list[HSP]:
    """Get a list of HSP objects that lack an alignment

    Args:
        gbks (list[GBK]): list of GBKs tos earch for missing gene hsp aligments

    Returns:
        list[HSP]: List of HSPs with missing alignments
    """
    hsps_to_align = []
    for gbk in gbks:
        for gene in gbk.genes:
            for hsp in gene.hsps:
                if hsp.alignment is None:
                    hsps_to_align.append(hsp)
    return hsps_to_align
