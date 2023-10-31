"""Contains enums relating to genbank files"""

from enum import Enum


class RECORD_TYPE(Enum):
    REGION = "region"
    CANDIDATE_CLUSTER = "cand_cluster"
    PROTO_CLUSTER = "proto_cluster"
    PROTO_CORE = "proto_core"
