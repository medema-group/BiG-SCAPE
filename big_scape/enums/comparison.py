"""Enum for comparison types."""

# from python
from enum import Enum


class ALIGNMENT_MODE(Enum):
    GLOBAL = "global"
    GLOCAL = "glocal"
    AUTO = "auto"


class COMPARISON_MODE(Enum):
    CDS = "cds"
    DOMAIN = "domain"


class LCS_MODE(Enum):
    REGION = "region"
    PROTOCLUSTER = "protocluster"
