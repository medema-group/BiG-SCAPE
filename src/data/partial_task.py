"""Contains code that scans the database for existing data and determines at which stage
of executing BiG-SCAPE the application last exited

Also contains functions to determine subsets of tasks that need to be done
"""

# from python
from __future__ import annotations
from enum import Enum
from typing import TYPE_CHECKING

# from other modules
from src.data import DB

# from circular imports
if TYPE_CHECKING:
    from src.genbank.gbk import GBK


class TASK(Enum):
    LOAD_GBKS = 0
    HMM_SCAN = 1
    HMM_ALIGN = 2
    COMPARISON = 3
    NOTHING_TO_DO = 4


def find_minimum_task(gbks: list[GBK]):
    """Finds the earliest task to start at. if new data was added, this will always
    be the load_gbks task. otherwise, it tries to find the latest task with unfinished
    business
    """
    input_data_state = get_input_data_state(gbks)

    if input_data_state.value < INPUT_TASK.SAME_DATA.value:
        return TASK.LOAD_GBKS

    hmm_data_state = get_hmm_data_state(gbks)

    if hmm_data_state.value < HMM_TASK.NEED_ALIGN.value:
        return TASK.HMM_SCAN

    if hmm_data_state.value < HMM_TASK.ALL_ALIGNED.value:
        return TASK.HMM_ALIGN

    comparison_data_state = get_comparison_data_state(gbks)

    if comparison_data_state.value < COMPARISON_TASK.ALL_DONE.value:
        return TASK.COMPARISON

    return TASK.NOTHING_TO_DO


class INPUT_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new data incoming, all data in database accounted for
    PARTIAL_DATA = 2  # some data that is in DB is not in input set
    MIXED_DATA = 3  # some data is new, some data is old
    SAME_DATA = 4  # data is the same


def get_input_data_state(gbks: list[GBK]) -> INPUT_TASK:
    """Returns the status of input data (gbks and regions) in the in-memory database"""
    distance_count = DB.get_table_row_count("gbk")

    if distance_count == 0:
        return INPUT_TASK.NO_DATA

    gbk_table = DB.metadata.tables["gbk"]

    # get set of gbks in database
    db_gbk_rows = DB.execute(gbk_table.select()).all()
    db_gbk_paths = set([db_gbk_row[1] for db_gbk_row in db_gbk_rows])
    input_gbk_paths = set([str(gbk.path) for gbk in gbks])

    if db_gbk_paths == input_gbk_paths:
        return INPUT_TASK.SAME_DATA

    sym_dif = db_gbk_paths.symmetric_difference(input_gbk_paths)

    # still same amount in db. new data
    if len(sym_dif) == len(db_gbk_paths):
        return INPUT_TASK.NEW_DATA

    # same amount in new data. there was more in db than in new data
    if len(sym_dif) == len(input_gbk_paths):
        return INPUT_TASK.PARTIAL_DATA

    # otherwise there is some new data, some old data is missing
    return INPUT_TASK.MIXED_DATA


class HMM_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEED_SCAN = 1  # new CDS need to be scanned
    NEED_ALIGN = 2  # new HSP need alignment
    ALL_ALIGNED = 3  # all HSP in database were aligned


def get_hmm_data_state(gbks: list[GBK]) -> HMM_TASK:
    """Retuns the state of data hmm processing in the in-memory database"""
    hsp_count = DB.get_table_row_count("hsp")

    if hsp_count == 0:
        return HMM_TASK.NO_DATA

    cds_count = DB.get_table_row_count("cds")
    scanned_cds_count = DB.get_table_row_count("scanned_cds")

    if cds_count > scanned_cds_count:
        return HMM_TASK.NEED_SCAN

    align_count = DB.get_table_row_count("hsp_alignment")

    if align_count < hsp_count:
        return HMM_TASK.NEED_ALIGN

    return HMM_TASK.ALL_ALIGNED


class COMPARISON_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new comparisons to be done
    ALL_DONE = 2  # all pairwise comparisons present


def get_comparison_data_state(gbks: list[GBK]) -> COMPARISON_TASK:
    """Retuns the state of pairwise comparison data in the in-memory database"""

    distance_count = DB.get_table_row_count("distance")

    if distance_count == 0:
        return COMPARISON_TASK.NO_DATA

    return None
