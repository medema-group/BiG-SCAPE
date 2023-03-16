"""Module containing code to load and store AntiSMASH regions"""

import logging
from typing import Dict, Optional

from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.genbank.candidate_cluster import CandidateCluster


class Region:
    """
    Class to describe a region within an Antismash GBK

    Attributes:
        number: int
        cand_clusters: Dict[int, CandidateCluster]
    """

    def __init__(self, number: int):
        self.number = number
        self.cand_clusters: Dict[int, Optional[CandidateCluster]] = {}

    def add_cand_cluster(self, cand_cluster: CandidateCluster):
        """Add a candidate cluster object to this region"""

        if cand_cluster.number not in self.cand_clusters:
            raise InvalidGBKRegionChildError()

        self.cand_clusters[cand_cluster.number] = cand_cluster

    @classmethod
    def parse(cls, feature: SeqFeature):
        """Creates a region object from a region feature in a GBK file"""
        if feature.type != "region":
            logging.error(
                "Feature is not of correct type! (expected: region, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "region_number" not in feature.qualifiers:
            logging.error("region_number qualifier not found in region feature!")
            raise InvalidGBKError()

        region_number = int(feature.qualifiers["region_number"][0])

        region = cls(region_number)

        if "candidate_cluster_numbers" not in feature.qualifiers:
            logging.error(
                "candidate_cluster_numbers qualifier not found in region feature!"
            )
            raise InvalidGBKError()

        for cand_cluster_number in feature.qualifiers["candidate_cluster_numbers"]:
            region.cand_clusters[int(cand_cluster_number)] = None

        return region
