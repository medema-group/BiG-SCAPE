"""Contains enums relating to genbank components"""

# from python
from enum import Enum


class COMPONENTS(Enum):
    GBK = "gbk"
    REGION = "region"
    PROTOCLUSTER = "protocluster"
    PROTO_CORE = "proto_core"
    CAND_CLUSTER = "cand_cluster"
    CDS = "CDS"
    CONCAT_CDS = "concat_cds"
