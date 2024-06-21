"""Enum for comparison types."""

# from python
from enum import Enum


class ALIGNMENT_MODE(Enum):
    GLOBAL = "global"
    GLOCAL = "glocal"
    LOCAL = "local"
    AUTO = "auto"


class COMPARISON_MODE(Enum):
    CDS = "cds"
    DOMAIN = "domain"


class LCS_MODE(Enum):
    REGION = "region"
    PROTOCLUSTER = "protocluster"


class CLASSIFY_MODE(Enum):
    CLASS = "class"
    CATEGORY = "category"
