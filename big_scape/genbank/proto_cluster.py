"""Module containing code to load and store AntiSMASH protoclusters"""

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
from big_scape.genbank.proto_core import ProtoCore


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import CandidateCluster
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files


class ProtoCluster(BGCRecord):
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        parent_gbk: GBK | None
        number: int
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
        categoty: str
        proto_core: dict[int, Optional[Protocore]]
        proto_core_cds_idx: set[int]

        _db_id: int | None
        _families: dict[float, int]
    """

    def __init__(
        self,
        parent_gbk: Optional[GBK],
        number: int,
        nt_start: int,
        nt_stop: int,
        contig_edge: Optional[bool],
        product: str,
        category: str,
        proto_core: dict[int, Optional[ProtoCore]],
    ):
        super().__init__(
            parent_gbk,
            number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )
        self.category: str = category
        self.proto_core: Dict[int, Optional[ProtoCore]] = proto_core
        self.proto_core_cds_idx: set[int] = set()

    def add_proto_core(self, proto_core: ProtoCore):
        """Add a proto_core object to this region

        Args:
            proto_core (ProtoCore): protocore object

        Raises:
            InvalidGBKRegionChildError: invalid child-parent relationship
        """

        if proto_core.number not in self.proto_core:
            raise InvalidGBKRegionChildError()

        proto_core.category = self.category

        self.proto_core[proto_core.number] = proto_core

        # protocores may exist in isolation
        if self.parent_gbk is None:
            return

        for idx, cds in enumerate(self.parent_gbk.genes):
            # skip to whatever is in this protocluster
            if cds.nt_start < self.nt_start:
                continue

            if cds.nt_stop > self.nt_stop:
                break

            if (
                cds.nt_start >= proto_core.nt_start
                and cds.nt_stop <= proto_core.nt_stop
            ):
                self.proto_core_cds_idx.add(idx)

    def save(self, parent_id: int, commit=True) -> None:
        """Stores this protocluster in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        super().save_record("protocluster", parent_id, commit)

    def save_all(self, parent_id: int) -> None:
        """Stores this protocluster and its children in the database. Does not
        commit immediately
        """
        self.save(parent_id, False)
        for protocore in self.proto_core.values():
            if protocore is None:
                continue

            if self._db_id is None:
                raise RuntimeError("Protocluster has no database id!")

            protocore.save(self._db_id, False)

    @classmethod
    def parse(
        cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None
    ) -> ProtoCluster:
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

        proto_core: dict[int, Optional[ProtoCore]] = {proto_cluster_number: None}

        category = ""
        if "category" in feature.qualifiers:
            category = feature.qualifiers["category"][0]

        nt_start, nt_stop, contig_edge, product = BGCRecord.parse_common(feature)

        return cls(
            parent_gbk,
            proto_cluster_number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
            category,
            proto_core,
        )

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

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

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
                record_table.c.product,
            )
            .where(record_table.c.record_type == "protocluster")
            .compile()
        )

        cursor_result = DB.execute(protocluster_select_query)

        protocluster_dict = {}

        for result in cursor_result.all():
            parent_candidate_cluster = candidate_cluster_dict[result.parent_id]
            parent_gbk = parent_candidate_cluster.parent_gbk

            new_proto_cluster = ProtoCluster(
                parent_gbk,
                result.record_number,
                result.nt_start,
                result.nt_stop,
                result.contig_edge,
                result.product,
                "",  # TODO: fix this
                {},
            )

            new_proto_cluster._db_id = result.id
            # add to parent CandidateCluster protocluster dict
            parent_candidate_cluster.proto_clusters[
                result.record_number
            ] = new_proto_cluster

            # add to dictionary
            protocluster_dict[result.id] = new_proto_cluster

        ProtoCore.load_all(protocluster_dict)
