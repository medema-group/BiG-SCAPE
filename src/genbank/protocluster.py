"""Module containing code to load and store AntiSMASH protoclusters"""

import logging

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError


class Protocluster:
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        number: int
        protocore: Protocore
    """

    def __init__(self, number: int, category: str):
        self.number = number
        self.category = category
        # self.protocore: Protocore = None

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

        if "category" not in feature.qualifiers:
            logging.error("category qualifier not found in protocluster feature!")
            raise InvalidGBKError()

        protocluster_category = feature.qualifiers["category"][0]

        protocluster = cls(protocluster_number, protocluster_category)

        return protocluster
