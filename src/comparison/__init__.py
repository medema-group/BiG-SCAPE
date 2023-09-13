"""Contains code to create and manage comparison objects"""
from .binning import (
    BGCPair,
    RecordPairGenerator,
    RecordPairGeneratorQueryRef,
    generate_mix,
)
from .comparable_region import ComparableRegion
from .legacy_workflow import create_bin_network_edges
from .legacy_workflow_alt import create_bin_network_edges_alt
from .legacy_bins import legacy_bin_generator, legacy_get_class

__all__ = [
    "BGCPair",
    "RecordPairGenerator",
    "RecordPairGeneratorQueryRef",
    "generate_mix",
    "ComparableRegion",
    "create_bin_network_edges",
    "create_bin_network_edges_alt",
    "legacy_bin_generator",
    "legacy_get_class",
]
