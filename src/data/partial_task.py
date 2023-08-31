"""Contains code that scans the database for existing data and determines at which stage
of executing BiG-SCAPE the application last exited

Also contains functions to determine subsets of tasks that need to be done 
"""

# from python
from enum import Enum

# from other modules
from src.data import DB


class TASK(Enum):
    LOAD_GBKS = 0
    HMM_SCAN = 1
    HMM_ALIGN = 2
    COMPARISON = 3
    NOTHING_TO_DO = 4


def find_minimum_task():
    """Finds the earliest task to start at. if new data was added, this will always
    be the load_gbks task. otherwise, it tries to find the latest task with unfinished
    business
    """
    input_data_state = get_input_data_state()

    if input_data_state.value < INPUT_TASK.SAME_DATA.value:
        return TASK.LOAD_GBKS

    hmm_data_state = get_hmm_data_state().value

    if hmm_data_state.value < HMM_TASK.ALL_SCANNED.value:
        return TASK.HMM_SCAN

    if hmm_data_state.value < HMM_TASK.ALL_ALIGNED.value:
        return TASK.HMM_SCAN

    comparison_data_state = get_comparison_data_state()

    if comparison_data_state.value < COMPARISON_TASK.ALL_DONE.value:
        return TASK.COMPARISON

    return TASK.NOTHING_TO_DO


class INPUT_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new data incoming, all data in database accounted for
    PARTIAL_DATA = 2  # some data that is in DB is not in input set
    SAME_DATA = 3  # data is the same


def get_input_data_state() -> INPUT_TASK:
    """Returns the status of input data (gbks and regions) in the in-memory database"""
    distance_count = DB.get_table_row_count("gbk")

    if distance_count == 0:
        return INPUT_TASK.NO_DATA


class HMM_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEED_SCAN = 0  # new CDS need to be scanned
    ALL_SCANNED = 1  # all CDS in database were scanned
    NEED_ALIGN = 2  # new HSP need alignment
    ALL_ALIGNED = 3  # all HSP in database were aligned


def get_hmm_data_state() -> HMM_TASK:
    """Retuns the state of data hmm processing in the in-memory database"""
    distance_count = DB.get_table_row_count("hsp")

    if distance_count == 0:
        return HMM_TASK.NO_DATA


class COMPARISON_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new comparisons to be done
    ALL_DONE = 2  # all pairwise comparisons present


def get_comparison_data_state() -> COMPARISON_TASK:
    """Retuns the state of pairwise comparison data in the in-memory database"""

    distance_count = DB.get_table_row_count("distance")

    if distance_count == 0:
        return COMPARISON_TASK.NO_DATA
