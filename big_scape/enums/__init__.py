"""Module containing code related to enums"""

from .input_parameters import INPUT_MODE
from .input_parameters import RUN_MODE
from .source_type import SOURCE_TYPE
from .partial_task import TASK, INPUT_TASK, HMM_TASK
from .comparison import (
    ALIGNMENT_MODE,
    EXTEND_STRATEGY,
    LCS_MODE,
    COMPARISON_MODE,
    CLASSIFY_MODE,
)
from .genbank import RECORD_TYPE, FEATURE_TYPE
from .components import COMPONENTS

__all__ = [
    "INPUT_MODE",
    "RUN_MODE",
    "SOURCE_TYPE",
    "TASK",
    "INPUT_TASK",
    "HMM_TASK",
    "ALIGNMENT_MODE",
    "EXTEND_STRATEGY",
    "LCS_MODE",
    "COMPARISON_MODE",
    "CLASSIFY_MODE",
    "RECORD_TYPE",
    "FEATURE_TYPE",
    "COMPONENTS",
]
