"""Module containing code to load and store AntiSMASH candidate clusters"""

# from python
from __future__ import annotations
import logging
from typing import Dict, Optional, TYPE_CHECKING

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from big_scape.data import DB
from big_scape.errors import InvalidGBKError, InvalidGBKRegionChildError

# from this module
from big_scape.genbank.bgc_record import BGCRecord
from big_scape.genbank.proto_cluster import ProtoCluster


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import Region  # imported earlier in file_input.load_files
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files


class CandidateCluster(BGCRecord):
    """
    Class to describe a candidate cluster within an Antismash GBK

    Attributes:
        parent_gbk: GBK | None
        number: int
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
        _db_id: int | None
        _families: dict[float, int]
        kind: str
        proto_clusters: Dict{number: int, ProtoCluster | None}
    """

    def __init__(
        self,
        parent_gbk: Optional[GBK],
        number: int,
        nt_start: int,
        nt_stop: int,
        contig_edge: Optional[bool],
        product: str,
        kind: str,
        proto_clusters: dict[int, Optional[ProtoCluster]],
    ):
        super().__init__(
            parent_gbk,
            number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )
        self.kind: str = kind
        self.proto_clusters: Dict[int, Optional[ProtoCluster]] = proto_clusters

    def add_proto_cluster(self, proto_cluster: ProtoCluster) -> None:
        """Add a protocluster object to this region

        Args:
            proto_cluster (ProtoCluster): antiSMASH protocluster

        Raises:
            InvalidGBKRegionChildError: invalid child-parent relationship
        """

        if proto_cluster.number not in self.proto_clusters:
            raise InvalidGBKRegionChildError()

        self.proto_clusters[proto_cluster.number] = proto_cluster

    def save(self, parent_id: int, commit=True) -> None:
        """Stores this candidate cluster in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        super().save_record("cand_cluster", parent_id, commit)

    def save_all(self, parent_id: int) -> None:
        """Stores this candidate cluster and its children in the database. Does not
        commit immediately
        """
        self.save(parent_id, False)
        for proto_cluster in self.proto_clusters.values():
            if proto_cluster is None:
                continue

            if self._db_id is None:
                raise RuntimeError("Candidate cluster has no database id!")

            proto_cluster.save_all(self._db_id)

    @classmethod
    def parse(
        cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None
    ) -> CandidateCluster:
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

        if "protoclusters" not in feature.qualifiers:
            logging.error("protoclusters qualifier not found in region feature!")
            raise InvalidGBKError()

        proto_clusters: dict[int, Optional[ProtoCluster]] = {}
        for proto_cluster_number in feature.qualifiers["protoclusters"]:
            proto_clusters[int(proto_cluster_number)] = None

        nt_start, nt_stop, contig_edge, product = BGCRecord.parse_common(feature)

        return cls(
            parent_gbk,
            cand_cluster_number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
            cand_cluster_kind,
            proto_clusters,
        )

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

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

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
                record_table.c.product,
            )
            .where(record_table.c.record_type == "cand_cluster")
            .where(record_table.c.parent_id.in_(region_dict.keys()))
            .compile()
        )

        cursor_result = DB.execute(candidate_cluster_select_query)

        candidate_cluster_dict = {}

        for result in cursor_result.all():
            parent_region = region_dict[result.parent_id]

            parent_gbk = parent_region.parent_gbk

            new_candidate_cluster = CandidateCluster(
                parent_gbk,
                result.record_number,
                result.nt_start,
                result.nt_stop,
                result.contig_edge,
                result.product,
                "",  # TODO: fix this
                {},
            )

            new_candidate_cluster._db_id = result.id

            # add to parent Region candidate cluster dict
            parent_region.cand_clusters[result.record_number] = new_candidate_cluster

            # add to dictionary
            candidate_cluster_dict[result.id] = new_candidate_cluster

        ProtoCluster.load_all(candidate_cluster_dict)
