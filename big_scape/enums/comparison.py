"""Enum for comparison types."""

# from python
from enum import Enum


class ALIGNMENT_MODE(Enum):
    GLOBAL = "global"
    GLOCAL = "glocal"
    LOCAL = "local"
    AUTO = "auto"


class EXTEND_STRATEGY(Enum):
    LEGACY = "legacy"
    GREEDY = "greedy"
    SIMPLE_MATCH = "simple_match"


class COMPARISON_MODE(Enum):
    CDS = "cds"
    DOMAIN = "domain"


class LCS_MODE(Enum):
    REGION = "region"
    PROTOCLUSTER = "protocluster"


class CLASSIFY_MODE(Enum):
    NONE = "none"
    CLASS = "class"
    CATEGORY = "category"
