"""Contains code to create and manage comparison objects"""
from .binning import BGCPair, BGCBin, generate_mix
from .comparable_region import ComparableRegion

__all__ = [
    "BGCPair",
    "BGCBin",
    "generate_mix",
    "ComparableRegion",
]
