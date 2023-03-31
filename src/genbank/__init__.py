"""Module containing code related to genbank file parsing and storage"""
from src.genbank.gbk import GBK, SOURCE_TYPE
from src.genbank.region import Region
from src.genbank.candidate_cluster import CandidateCluster
from src.genbank.proto_cluster import ProtoCluster
from src.genbank.proto_core import ProtoCore
from src.genbank.cds import CDS
from src.genbank.bgc_record import BGCRecord

__all__ = [
    "GBK",
    "Region",
    "CandidateCluster",
    "ProtoCluster",
    "ProtoCore",
    "CDS",
    "BGCRecord",
    "SOURCE_TYPE",
]
