"""Contains code to create and manage comparison objects"""
from .binning import (
    RecordPair,
    RecordPairGenerator,
    QueryToRefRecordPairGenerator,
    RefToRefRecordPairGenerator,
    MissingRecordPairGenerator,
    generate_mix,
)
from .comparable_region import ComparableRegion
from .legacy_workflow_alt import generate_edges
from .legacy_bins import legacy_bin_generator, legacy_get_class
from .utility import save_edge_to_db

from . import lcs

__all__ = [
    "RecordPair",
    "RecordPairGenerator",
    "QueryToRefRecordPairGenerator",
    "RefToRefRecordPairGenerator",
    "MissingRecordPairGenerator",
    "generate_mix",
    "ComparableRegion",
    "generate_edges",
    "legacy_bin_generator",
    "legacy_get_class",
    "save_edge_to_db",
    "lcs",
]
