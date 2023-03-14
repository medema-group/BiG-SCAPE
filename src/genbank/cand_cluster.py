"""Module containing code to load and store AntiSMASH candidate clusters"""

import logging
from typing import Dict

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError
from src.genbank.protocluster import Protocluster


class CandidateCluster:
    """
    Class to describe a candidate cluster within an Antismash GBK

    Attributes:
        number: int
        kind: str
        protoclusters: Dict[int, Protocluster]
    """

    def __init__(self, number: int, kind: str):
        self.number = number
        self.kind = kind
        self.protoclusters: Dict[int, Protocluster] = {}

    @classmethod
    def parse(cls, feature: SeqFeature):
        """Creates a cand_cluster object from a region feature in a GBK file"""
        if feature.type != "cand_cluster":
            logging.error(
                "Feature is not of correct type! (expected: cand_cluster, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "candidate_cluster_number" not in feature.qualifiers:
            logging.error(
                "candidate_cluster_number qualifier not found in cand_cluster feature!"
            )
            raise InvalidGBKError()

        cand_cluster_number = int(feature.qualifiers["candidate_cluster_number"][0])

        if "kind" not in feature.qualifiers:
            logging.error("kind qualifier not found in cand_cluster feature!")
            raise InvalidGBKError()

        cand_cluster_kind = feature.qualifiers["kind"][0]

        cand_cluster = cls(cand_cluster_number, cand_cluster_kind)

        return cand_cluster
