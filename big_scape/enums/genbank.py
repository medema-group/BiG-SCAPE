"""Contains enums relating to genbank files"""

from enum import Enum


class RECORD_TYPE(Enum):
    REGION = "region"
    CAND_CLUSTER = "cand_cluster"
    PROTO_CLUSTER = "protocluster"
    PROTO_CORE = "proto_core"


class FEATURE_TYPE(Enum):
    REGION = "region"
    PROTOCLUSTER = "protocluster"
    PROTO_CORE = "proto_core"
    CAND_CLUSTER = "cand_cluster"
    CDS = "CDS"
    CLUSTER = "cluster"
