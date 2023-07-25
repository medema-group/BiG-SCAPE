"""Module containing code to load and store AntiSMASH regions"""

# from python
from __future__ import annotations
import logging
from typing import Dict, Optional, TYPE_CHECKING

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB
from src.errors import InvalidGBKError, InvalidGBKRegionChildError

# from this module
from src.genbank.bgc_record import BGCRecord
from src.genbank.candidate_cluster import CandidateCluster
from src.genbank.cds import CDS

# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


class Region(BGCRecord):
    """
    Class to describe a region within an Antismash GBK

    Attributes:
        number: int
        cand_clusters: Dict{number: int, CandidateCluster}
    """

    def __init__(self, number: int):
        super().__init__(number)
        self.cand_clusters: Dict[int, Optional[CandidateCluster]] = {}

    def add_cand_cluster(self, cand_cluster: CandidateCluster) -> None:
        """Add a candidate cluster object to this region

        Args:
            cand_cluster (CandidateCluster): candidate cluster object

        Raises:
            InvalidGBKRegionChildError: Invalid gbk region child
        """

        if cand_cluster.number not in self.cand_clusters:
            raise InvalidGBKRegionChildError()

        self.cand_clusters[cand_cluster.number] = cand_cluster

    def save(self, commit=True) -> None:
        """Stores this region in the database

        Args:
            commit: commit immediately after executing the insert query"""
        super().save_record("region", commit)

    def save_all(self) -> None:
        """Stores this Region and its children in the database. Does not commit immediately"""
        self.save(False)
        for candidate_cluster in self.cand_clusters.values():
            if candidate_cluster is None:
                continue

            candidate_cluster.save_all()

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None) -> Region:
        """Creates a region object from a region feature in a GBK file

        Args:
            feature (SeqFeature): region(as5+) or cluster (as4) GBK feature

        Raises:
            InvalidGBKError: Invalid or missing fields

        Returns:
            Region: region object
        """
        if feature.type != "region" and feature.type != "cluster":
            logging.error(
                "Feature is not of correct type! (expected: region or cluster, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        # AS5 and up gbks have region features, as well as candidate clusters and
        # children classes (protocluster, protocore)
        if feature.type == "region":
            if "region_number" not in feature.qualifiers:
                logging.error("region number qualifier not found in region feature!")
                raise InvalidGBKError()

            region_number = int(feature.qualifiers["region_number"][0])

            region = cls(region_number)

            region.parse_bgc_record(feature, parent_gbk=parent_gbk)

            if "candidate_cluster_numbers" not in feature.qualifiers:
                logging.error(
                    "candidate_cluster_numbers qualifier not found in region feature!"
                )
                raise InvalidGBKError()

            for cand_cluster_number in feature.qualifiers["candidate_cluster_numbers"]:
                region.cand_clusters[int(cand_cluster_number)] = None

            return region

        # AS4 gbks have cluster features instead of region, and no children features
        # we artifically input the info in the cluster feature into the Region object
        if feature.type == "cluster":
            if (
                "note" not in feature.qualifiers
                or "Cluster number" not in feature.qualifiers["note"][0]
            ):
                logging.error("cluster number qualifier not found in cluster feature!")
                raise InvalidGBKError()

            cluster_note_number = feature.qualifiers["note"][0]
            cluster_number = int(cluster_note_number.split(": ")[1])
            region = cls(cluster_number)

            region.parse_bgc_record(feature, parent_gbk=parent_gbk)
            return region

        else:
            raise ValueError("Could not parse region feature")

    def get_cds_with_domains(self, return_all=True, reverse=False) -> list[CDS]:
        return super().get_cds_with_domains(return_all, reverse)

    def get_attr_dict(self) -> dict[str, object]:
        """Gets a dictionary of attributes, useful for adding to network nodes later"""
        return super().get_attr_dict()

    def __repr__(self):
        return f"{self.parent_gbk} Region {self.number} {self.nt_start}-{self.nt_stop} "

    @staticmethod
    def load_all(gbk_dict: dict[int, GBK]) -> None:
        """Load all Region objects from the database

        This function populates the region objects in the GBKs provided in the input
        gbk_dict

        Args:
            gbk_dict (dict[int, GBK]): Dictionary of GBK objects with database ids as
            keys. Used for reassembling the hierarchy
        """
        record_table = DB.metadata.tables["bgc_record"]

        region_select_query = (
            record_table.select()
            .add_columns(
                record_table.c.id,
                record_table.c.gbk_id,
                record_table.c.parent_id,
                record_table.c.record_type,
                record_table.c.contig_edge,
                record_table.c.nt_start,
                record_table.c.nt_stop,
                record_table.c.product,
            )
            .where(record_table.c.record_type == "region")
            .compile()
        )

        cursor_result = DB.execute(region_select_query)

        region_dict = {}

        for result in cursor_result.all():
            new_region = Region(1)

            new_region._db_id = result.id

            new_region.parent_gbk = gbk_dict[result.gbk_id]
            new_region.nt_start = result.nt_start
            new_region.nt_stop = result.nt_stop
            new_region.contig_edge = result.contig_edge
            new_region.product = result.product

            # add to parent GBK
            gbk_dict[result.gbk_id].region = new_region

            # add to dictionary
            region_dict[result.id] = new_region

        CandidateCluster.load_all(region_dict)
