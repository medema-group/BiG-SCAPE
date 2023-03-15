"""Module containing code to load and store AntiSMASH protoclusters"""

import logging
from typing import Dict, Optional

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.genbank.proto_core import Protocore


class Protocluster:
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
        self.protocore: Dict[int, Optional[Protocore]] = {}

    def add_proto_core(self, proto_core: Protocore):
        """Add a proto_core object to this region"""

        if proto_core.number not in self.protocore:
            raise InvalidGBKRegionChildError()

        self.protocore[proto_core.number] = proto_core

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

        protocluster_number = int(feature.qualifiers["protocluster_number"][0])

        protocluster = cls(protocluster_number)
        protocluster.protocore[protocluster_number] = None

        if "category" not in feature.qualifiers:
            logging.error("category qualifier not found in protocluster feature!")
            raise InvalidGBKError()

        protocluster_category = feature.qualifiers["category"][0]

        protocluster.category = protocluster_category

        return protocluster
