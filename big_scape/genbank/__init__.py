"""Module containing code related to genbank file parsing and storage"""
from .gbk import GBK
from .region import Region
from .candidate_cluster import CandidateCluster
from .proto_cluster import ProtoCluster
from .proto_core import ProtoCore, MergedProtoCore
from .cds import (
    CDS,
    check_translation,
    get_translation,
    translate,
    trim_fuzzy,
)
from .bgc_record import BGCRecord

__all__ = [
    "GBK",
    "Region",
    "CandidateCluster",
    "ProtoCluster",
    "ProtoCore",
    "MergedProtoCore",
    "CDS",
    "check_translation",
    "get_translation",
    "translate",
    "trim_fuzzy",
    "BGCRecord",
]
