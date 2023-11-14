"""Contains enums for partial task detection"""

from enum import Enum


class TASK(Enum):
    SAVE_GBKS = 0
    HMM_SCAN = 1
    HMM_ALIGN = 2
    COMPARISON = 3
    NOTHING_TO_DO = 4


class INPUT_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new data incoming, all data in database accounted for
    PARTIAL_DATA = 2  # some data that is in DB is not in input set
    MIXED_DATA = 3  # some data is new, some data is old
    SAME_DATA = 4  # data is the same


class HMM_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEED_SCAN = 1  # new CDS need to be scanned
    NEED_ALIGN = 2  # new HSP need alignment
    ALL_ALIGNED = 3  # all HSP in database were aligned


class COMPARISON_TASK(Enum):
    NO_DATA = 0  # nothing done yet
    NEW_DATA = 1  # new comparisons to be done
    ALL_DONE = 2  # all pairwise comparisons present
