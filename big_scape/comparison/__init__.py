"""Contains code to create and manage comparison objects"""

from .binning import (
    RecordPairGenerator,
    MissingRecordPairGenerator,
    QueryRecordPairGenerator,
    QueryMissingRecordPairGenerator,
    ConnectedComponentPairGenerator,
    generate_mix_bin,
    legacy_bin_generator,
    legacy_get_class,
    as_class_bin_generator,
    get_legacy_weights_from_category,
    get_record_category,
)
from .comparable_region import ComparableRegion
from .workflow import generate_edges
from .utility import (
    save_edge_to_db,
    save_edges_to_db,
    get_edge_param_id,
    edge_params_insert,
    edge_params_query,
    get_edge_weight,
)

from . import lcs
from . import extend

__all__ = [
    "record_pair",
    "RecordPairGenerator",
    "MissingRecordPairGenerator",
    "QueryRecordPairGenerator",
    "QueryMissingRecordPairGenerator",
    "ConnectedComponentPairGenerator",
    "generate_mix_bin",
    "ComparableRegion",
    "generate_edges",
    "legacy_bin_generator",
    "legacy_get_class",
    "as_class_bin_generator",
    "get_legacy_weights_from_category",
    "get_record_category",
    "save_edge_to_db",
    "save_edges_to_db",
    "lcs",
    "extend",
    "get_edge_param_id",
    "edge_params_insert",
    "edge_params_query",
    "get_edge_weight",
]
