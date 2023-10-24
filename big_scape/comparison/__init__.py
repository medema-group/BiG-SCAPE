"""Contains code to create and manage comparison objects"""
from .binning import (
    RecordPair,
    RecordPairGenerator,
    QueryToRefRecordPairGenerator,
    RefToRefRecordPairGenerator,
    MissingRecordPairGenerator,
    generate_mix,
    legacy_bin_generator,
    legacy_get_class,
    as_class_bin_generator,
    get_weight_category,
    get_record_category,
)
from .comparable_region import ComparableRegion
from .legacy_workflow_alt import generate_edges
from .utility import save_edge_to_db, save_edges_to_db

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
    "as_class_bin_generator",
    "get_weight_category",
    "get_record_category",
    "save_edge_to_db",
    "lcs",
    "save_edges_to_db",
]
