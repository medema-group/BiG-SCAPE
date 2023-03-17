"""Module containing code to load and store AntiSMASH protoclusters"""

# from python
import logging
from typing import Dict, Optional

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError, InvalidGBKRegionChildError

# from this module
from src.genbank.bgc_record import BGCRecord
from src.genbank.proto_core import ProtoCore


class ProtoCluster(BGCRecord):
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        number: int
        categoty: str
        protocore: Protocore
    """

    def __init__(self, number: int):
        super().__init__()
        self.number = number
        self.category: str = ""
        self.proto_core: Dict[int, Optional[ProtoCore]] = {}

    def add_proto_core(self, proto_core: ProtoCore):
        """Add a proto_core object to this region"""

        if proto_core.number not in self.proto_core:
            raise InvalidGBKRegionChildError()

        self.proto_core[proto_core.number] = proto_core

    def save(self, commit=True):
        """Stores this protocluster in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        return super().save("protocluster", commit)

    @classmethod
    def parse(cls, feature: SeqFeature):
        """Creates a Protocluster object from a region feature in a GBK file"""
        if feature.type != "protocluster":
            logging.error(
                "Feature is not of correct type! (expected: protocluster, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "protocluster_number" not in feature.qualifiers:
            logging.error(
                "protocluster_number qualifier not found in protocluster feature!"
            )
            raise InvalidGBKError()

        proto_cluster_number = int(feature.qualifiers["protocluster_number"][0])

        proto_cluster = cls(proto_cluster_number)
        proto_cluster.parse_location(feature)
        proto_cluster.proto_core[proto_cluster_number] = None

        if "category" not in feature.qualifiers:
            logging.error("category qualifier not found in protocluster feature!")
            raise InvalidGBKError()

        protocluster_category = feature.qualifiers["category"][0]

        proto_cluster.category = protocluster_category

        return proto_cluster
