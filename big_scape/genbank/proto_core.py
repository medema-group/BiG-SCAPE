"""Module containing code to load and store AntiSMASH protoclusters"""

# from python
from __future__ import annotations
import logging
from typing import Optional, TYPE_CHECKING

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from big_scape.data import DB
from big_scape.errors import InvalidGBKError

# from this module
from big_scape.genbank.bgc_record import BGCRecord


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import ProtoCluster
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files


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
        category: Optional[str] = None,
    ):
        super().__init__(
            parent_gbk,
            number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )

        self.category: Optional[str] = category

    def save(self, parent_id: int, commit=True) -> None:
        """Stores this protocore in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        super().save_record("proto_core", parent_id, commit)

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
            protocluster_dict (dict[int, GBK]): Dictionary of protocluster objects with database ids
            as keys. Used for parenting
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

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
                record_table.c.category,
                record_table.c.merged,
            )
            .where(record_table.c.record_type == "proto_core")
            .where(record_table.c.parent_id.in_(protocluster_dict.keys()))
            .compile()
        )

        cursor_result = DB.execute(region_select_query)

        for result in cursor_result.all():
            parent_proto_cluster = protocluster_dict[result.parent_id]
            parent_gbk = parent_proto_cluster.parent_gbk

            merged = result.merged
            # if merged = True load as MergedProtoCore and load merged number as string
            # if false load as normal protocore
            if merged:
                new_proto_core: Optional[ProtoCore] = MergedProtoCore(
                    parent_gbk,
                    result.record_number,
                    result.nt_start,
                    result.nt_stop,
                    result.contig_edge,
                    result.product,
                    result.category,
                )

            else:
                new_proto_core = ProtoCore(
                    parent_gbk,
                    result.record_number,
                    result.nt_start,
                    result.nt_stop,
                    result.contig_edge,
                    result.product,
                    result.category,
                )

            if new_proto_core is not None:
                new_proto_core._db_id = result.id

            # add to parent ProtoCluster protocore dict
            parent_proto_cluster.proto_core[result.record_number] = new_proto_core


class MergedProtoCore(ProtoCore):
    """Class to described a merged protocore within an Antismash GBK

    Args:
        ProtoCore (_type_): _description_
    """

    def __init__(
        self,
        parent_gbk: Optional[GBK],
        merged_number: str,
        nt_start: int,
        nt_stop: int,
        contig_edge: Optional[bool],
        product: str,
        category: Optional[str] = None,
    ):
        super().__init__(
            parent_gbk,
            0,
            nt_start,
            nt_stop,
            contig_edge,
            product,
            category,
        )

        self.merged_number = merged_number
        self.merged = True

    @staticmethod
    def merge(protocore_a, protocore_b) -> MergedProtoCore:
        """Merges two protocores into a single merged protocore

        Args:
            protocore_a (_type_): _description_
            protocore_b (_type_): _description_

        Returns:
            MergedProtoCore: _description_
        """

        if protocore_a.parent_gbk != protocore_b.parent_gbk:
            raise ValueError("Cannot merge protocores from different GBKs")

        parent_gbk = protocore_a.parent_gbk
        merged_number = f"{protocore_a.number}_{protocore_b.number}"

        if protocore_a.category != protocore_b.category:
            category = f"{protocore_a.category}.{protocore_b.category}"
        else:
            category = protocore_a.category

        if protocore_a.product != protocore_b.product:
            product = f"{protocore_a.product}.{protocore_b.product}"
        else:
            product = protocore_a.product

        contig_edge = protocore_a.contig_edge or protocore_b.contig_edge

        nt_start = min(protocore_a.nt_start, protocore_b.nt_start)
        nt_stop = max(protocore_a.nt_stop, protocore_b.nt_stop)

        merged_proto_core = MergedProtoCore(
            parent_gbk,
            merged_number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
            category,
        )

        return merged_proto_core
