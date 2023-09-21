"""Contains code to create and manage comparison objects"""
from .binning import (
    RecordPair,
    RecordPairGenerator,
    RecordPairGeneratorQueryRef,
    RecordPairGeneratorConRefSinRef,
    PartialRecordPairGenerator,
    generate_mix,
)
from .comparable_region import ComparableRegion
from .legacy_workflow import create_bin_network_edges
from .legacy_workflow_alt import create_bin_network_edges_alt, generate_edges
from .legacy_bins import legacy_bin_generator, legacy_get_class
from .utility import save_edge_to_db

__all__ = [
    "RecordPair",
    "RecordPairGenerator",
    "RecordPairGeneratorQueryRef",
    "RecordPairGeneratorConRefSinRef",
    "PartialRecordPairGenerator",
    "generate_mix",
    "ComparableRegion",
    "create_bin_network_edges",
    "create_bin_network_edges_alt",
    "generate_edges",
    "legacy_bin_generator",
    "legacy_get_class",
    "save_edge_to_db",
]
