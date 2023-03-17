"""Module containing code to load and store AntiSMASH regions"""

# from python
import logging
from typing import Dict, Optional

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError, InvalidGBKRegionChildError

# from this module
from src.genbank.bgc_record import BGCRecord
from src.genbank.candidate_cluster import CandidateCluster


class Region(BGCRecord):
    """
    Class to describe a region within an Antismash GBK

    Attributes:
        number: int
        cand_clusters: Dict[int, CandidateCluster]
    """

    def __init__(self, number: int):
        super().__init__()
        self.number = number
        self.cand_clusters: Dict[int, Optional[CandidateCluster]] = {}

    def add_cand_cluster(self, cand_cluster: CandidateCluster):
        """Add a candidate cluster object to this region"""

        if cand_cluster.number not in self.cand_clusters:
            raise InvalidGBKRegionChildError()

        self.cand_clusters[cand_cluster.number] = cand_cluster

    def save(self, commit=True):
        """Stores this region in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        return super().save("region", commit)

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

        region.parse_location(feature)

        if "candidate_cluster_numbers" not in feature.qualifiers:
            logging.error(
                "candidate_cluster_numbers qualifier not found in region feature!"
            )
            raise InvalidGBKError()

        for cand_cluster_number in feature.qualifiers["candidate_cluster_numbers"]:
            region.cand_clusters[int(cand_cluster_number)] = None

        return region
