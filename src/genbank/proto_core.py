"""Module containing code to load and store AntiSMASH protoclusters"""

# from python
import logging

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError


class ProtoCore:
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        number: int
    """

    def __init__(self, number: int):
        self.number = number

    @classmethod
    def parse(cls, feature: SeqFeature):
        """Creates a Protocore object from a region feature in a GBK file"""
        if feature.type != "proto_core":
            logging.error(
                "Feature is not of correct type! (expected: proto_core, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "protocluster_number" not in feature.qualifiers:
            logging.error(
                "protocluster_number qualifier not found in proto_core feature!"
            )
            raise InvalidGBKError()

        proto_core_number = int(feature.qualifiers["protocluster_number"][0])

        proto_core = cls(proto_core_number)

        return proto_core
