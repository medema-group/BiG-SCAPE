"""Module containing code to load and store AntiSMASH candidate clusters"""

import logging
from typing import Dict, Optional

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.genbank.proto_cluster import Protocluster


class CandidateCluster:
    """
    Class to describe a candidate cluster within an Antismash GBK

    Attributes:
        number: int
        kind: str
        proto_clusters: Dict[int, Protocluster]
    """

    def __init__(self, number: int):
        self.number = number
        self.kind: str = ""
        self.proto_clusters: Dict[int, Optional[Protocluster]] = {}

    def add_proto_cluster(self, proto_cluster: Protocluster):
        """Add a protocluster object to this region"""

        if proto_cluster.number not in self.proto_clusters:
            raise InvalidGBKRegionChildError()

        self.proto_clusters[proto_cluster.number] = proto_cluster

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

        cand_cluster = cls(cand_cluster_number)
        cand_cluster.kind = cand_cluster_kind

        if "protoclusters" not in feature.qualifiers:
            logging.error("protoclusters qualifier not found in region feature!")
            raise InvalidGBKError()

        for proto_cluster_number in feature.qualifiers["protoclusters"]:
            cand_cluster.proto_clusters[int(proto_cluster_number)] = None

        return cand_cluster
