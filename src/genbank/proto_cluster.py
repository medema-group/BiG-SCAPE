"""Module containing code to load and store AntiSMASH protoclusters"""

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
from src.genbank.proto_core import ProtoCore


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import CandidateCluster
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


class ProtoCluster(BGCRecord):
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
        number: int
        categoty: str
        protocore: Protocore
    """

    def __init__(self, number: int):
        super().__init__(number)
        self.category: str = ""
        self.proto_core: Dict[int, Optional[ProtoCore]] = {}

    def add_proto_core(self, proto_core: ProtoCore):
        """Add a proto_core object to this region

        Args:
            proto_core (ProtoCore): protocore object

        Raises:
            InvalidGBKRegionChildError: invalid child-parent relationship
        """

        if proto_core.number not in self.proto_core:
            raise InvalidGBKRegionChildError()

        self.proto_core[proto_core.number] = proto_core

    def save(self, commit=True):
        """Stores this protocluster in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        return super().save("protocluster", commit)

    def save_all(self):
        """Stores this protocluster and its children in the database. Does not
        commit immediately
        """
        self.save(False)
        for protocore in self.proto_core.values():
            protocore.save(False)

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None):
        """Creates a Protocluster object from a region feature in a GBK file

        Args:
            feature (SeqFeature): protocluster genbank feature

        Raises:
            InvalidGBKError: invalid or missing fields

        Returns:
            ProtoCluster: protocluster object
        """
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
        proto_cluster.parse_bgc_record(feature, parent_gbk=parent_gbk)
        proto_cluster.proto_core[proto_cluster_number] = None

        if "category" in feature.qualifiers:
            proto_cluster.category = feature.qualifiers["category"][0]

        return proto_cluster

    def __repr__(self) -> str:
        return f"{self.parent_gbk} ProtoCluster {self.number} {self.nt_start}-{self.nt_stop} "

    @staticmethod
    def load_all(candidate_cluster_dict: dict[int, CandidateCluster]):
        """Load all ProtoCluster objects from the database

        This function populates the CandidateCluster objects in the GBKs provided in the
        input candidate_cluster_dict

        Args:
            candidate_cluster_dict (dict[int, GBK]): Dictionary of CandidateCluster
            objects with database ids as keys. Used for reassembling the hierarchy
        """
        record_table = DB.metadata.tables["bgc_record"]

        protocluster_select_query = (
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
            .where(record_table.c.record_type == "protocluster")
            .compile()
        )

        cursor_result = DB.execute(protocluster_select_query)

        protocluster_dict = {}

        for result in cursor_result.all():
            new_proto_cluster = ProtoCluster(result.record_number)
            new_proto_cluster.nt_start = result.nt_start
            new_proto_cluster.nt_stop = result.nt_stop
            new_proto_cluster.contig_edge = result.contig_edge

            # add to parent CandidateCluster protocluster dict
            candidate_cluster_dict[result.parent_id].proto_clusters[
                result.record_number
            ] = new_proto_cluster

            # add to dictionary
            protocluster_dict[result.id] = new_proto_cluster

        ProtoCore.load_all(protocluster_dict)
