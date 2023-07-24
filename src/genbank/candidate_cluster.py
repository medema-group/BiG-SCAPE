"""Module containing code to load and store AntiSMASH candidate clusters"""

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
from src.genbank.proto_cluster import ProtoCluster


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import Region  # imported earlier in src.file_input.load_files
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


class CandidateCluster(BGCRecord):
    """
    Class to describe a candidate cluster within an Antismash GBK

    Attributes:
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
        number: int
        kind: str
        proto_clusters: Dict{number: int, ProtoCluster}
    """

    def __init__(self, number: int):
        super().__init__(number)
        self.kind: str = ""
        self.proto_clusters: Dict[int, Optional[ProtoCluster]] = {}

    def add_proto_cluster(self, proto_cluster: ProtoCluster):
        """Add a protocluster object to this region

        Args:
            proto_cluster (ProtoCluster): antiSMASH protocluster

        Raises:
            InvalidGBKRegionChildError: invalid child-parent relationship
        """

        if proto_cluster.number not in self.proto_clusters:
            raise InvalidGBKRegionChildError()

        self.proto_clusters[proto_cluster.number] = proto_cluster

    def save(self, commit=True):
        """Stores this candidate cluster in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        return super().save("cand_cluster", commit)

    def save_all(self):
        """Stores this candidate cluster and its children in the database. Does not
        commit immediately
        """
        self.save(False)
        for proto_cluster in self.proto_clusters.values():
            proto_cluster.save_all()

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None):
        """_summary_Creates a cand_cluster object from a region feature in a GBK file

        Args:
            feature (SeqFeature): cand_cluster GBK feature

        Raises:
            InvalidGBKError: invalid or missing fields

        Returns:
            CandidateCluster: Candidate cluster object
        """
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
        cand_cluster.parse_bgc_record(feature, parent_gbk=parent_gbk)
        cand_cluster.kind = cand_cluster_kind

        if "protoclusters" not in feature.qualifiers:
            logging.error("protoclusters qualifier not found in region feature!")
            raise InvalidGBKError()

        for proto_cluster_number in feature.qualifiers["protoclusters"]:
            cand_cluster.proto_clusters[int(proto_cluster_number)] = None

        return cand_cluster

    def __repr__(self) -> str:
        return f"{self.parent_gbk} Candidate cluster {self.number} {self.nt_start}-{self.nt_stop} "

    @staticmethod
    def load_all(region_dict: dict[int, Region]):
        """Load all CandidateCluster objects from the database

        This function populates the CandidateCluster lists in the Regions provided in
        the input region_dict

        Args:
            region_dict (dict[int, Region]): Dictionary of Region objects with database
            ids as keys. Used for reassembling the hierarchy
        """
        record_table = DB.metadata.tables["bgc_record"]

        candidate_cluster_select_query = (
            record_table.select()
            .add_columns(
                record_table.c.id,
                record_table.c.record_number,
                record_table.c.parent_id,
                record_table.c.record_type,
                record_table.c.contig_edge,
                record_table.c.nt_start,
                record_table.c.nt_stop,
            )
            .where(record_table.c.record_type == "cand_cluster")
            .compile()
        )

        cursor_result = DB.execute(candidate_cluster_select_query)

        candidate_cluster_dict = {}

        for result in cursor_result.all():
            new_candidate_cluster = CandidateCluster(result.record_number)

            new_candidate_cluster._db_id = result.id

            new_candidate_cluster.nt_start = result.nt_start
            new_candidate_cluster.nt_stop = result.nt_stop
            new_candidate_cluster.contig_edge = result.contig_edge

            # add to parent Region candidate cluster dict
            region_dict[result.parent_id].cand_clusters[
                result.record_number
            ] = new_candidate_cluster

            # add to dictionary
            candidate_cluster_dict[result.id] = new_candidate_cluster

        ProtoCluster.load_all(candidate_cluster_dict)
