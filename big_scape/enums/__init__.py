"""Module containing code related to enums"""

from .input_parameters import INPUT_MODE
from .source_type import SOURCE_TYPE
from .partial_task import TASK, INPUT_TASK, HMM_TASK, COMPARISON_TASK
from .comparison import (
    ALIGNMENT_MODE,
    EXTEND_STRATEGY,
    LCS_MODE,
    COMPARISON_MODE,
    CLASSIFY_MODE,
)
from .genbank import RECORD_TYPE

__all__ = [
    "INPUT_MODE",
    "SOURCE_TYPE",
    "TASK",
    "INPUT_TASK",
    "HMM_TASK",
    "COMPARISON_TASK",
    "ALIGNMENT_MODE",
    "EXTEND_STRATEGY",
    "LCS_MODE",
    "COMPARISON_MODE",
    "CLASSIFY_MODE",
    "RECORD_TYPE",
]
