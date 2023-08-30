"""Module containing code to load and store AntiSMASH protoclusters"""

# from python
from __future__ import annotations
import logging
from typing import Optional, TYPE_CHECKING

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB
from src.errors import InvalidGBKError

# from this module
from src.genbank.bgc_record import BGCRecord


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import ProtoCluster
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


class ProtoCore(BGCRecord):
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        parent_gbk: GBK | None
        number: int
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
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
    ):
        super().__init__(
            parent_gbk,
            number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )

    def save(self, commit=True) -> None:
        """Stores this protocore in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        super().save_record("proto_core", commit)

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None) -> ProtoCore:
        """Creates a Protocore object from a region feature in a GBK file

        Args:
            feature (SeqFeature): proto_core antiSMASH genbank feature

        Raises:
            InvalidGBKError: invalid or missing field

        Returns:
            ProtoCore: protocore object
        """
        if feature.type != "proto_core":
            logging.error(
                "Feature is not of correct type! (expected: proto_core, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "protocluster_number" not in feature.qualifiers:
            logging.error(
                "protocluster_number qualifier not found in proto_core feature!"
            )
            raise InvalidGBKError()

        proto_core_number = int(feature.qualifiers["protocluster_number"][0])

        nt_start, nt_stop, contig_edge, product = BGCRecord.parse_common(feature)

        return cls(
            parent_gbk, proto_core_number, nt_start, nt_stop, contig_edge, product
        )

    def __repr__(self) -> str:
        return (
            f"{self.parent_gbk} ProtoCore {self.number} {self.nt_start}-{self.nt_stop} "
        )

    @staticmethod
    def load_all(protocluster_dict: dict[int, ProtoCluster]):
        """Load all ProtoCore objects from the database

        This function populates the region objects in the GBKs provided in the input
        gbk_dict

        Args:
            region_dict (dict[int, GBK]): Dictionary of Region objects with database ids
            as keys. Used for parenting
        """
        record_table = DB.metadata.tables["bgc_record"]

        region_select_query = (
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

        cursor_result = DB.execute(region_select_query)

        for result in cursor_result.all():
            parent_proto_cluster = protocluster_dict[result.parent_id]
            parent_gbk = parent_proto_cluster.parent_gbk

            new_proto_core = ProtoCore(
                parent_gbk,
                result.record_number,
                result.nt_start,
                result.nt_stop,
                result.contig_edge,
                result.product,
            )

            new_proto_core._db_id = result.id

            # add to parent ProtoCluster protocore dict
            parent_proto_cluster.proto_core[result.record_number] = new_proto_core
