"""Module containing code to load and store AntiSMASH protoclusters"""

import logging
from typing import Dict, Optional

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.genbank.proto_core import ProtoCore


class ProtoCluster:
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        number: int
        categoty: str
        protocore: Protocore
    """

    def __init__(self, number: int):
        self.number = number
        self.category: str = ""
        self.proto_core: Dict[int, Optional[ProtoCore]] = {}

    def add_proto_core(self, proto_core: ProtoCore):
        """Add a proto_core object to this region"""

        if proto_core.number not in self.proto_core:
            raise InvalidGBKRegionChildError()

        self.proto_core[proto_core.number] = proto_core

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
        proto_cluster.proto_core[proto_cluster_number] = None

        if "category" not in feature.qualifiers:
            logging.error("category qualifier not found in protocluster feature!")
            raise InvalidGBKError()

        protocluster_category = feature.qualifiers["category"][0]

        proto_cluster.category = protocluster_category

        return proto_cluster
