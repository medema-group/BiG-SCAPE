"""Contains a class for a BGC record, a common set of qualifiers/attributes amongst the
AntiSMASh genbank records
"""

# from python
from __future__ import annotations
from typing import Any, Optional
import logging

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB
from src.errors import InvalidGBKError

# from this module


class BGCRecord:
    """Describes a common set of qualifiers/attributes among all ඞ AntiSMASH genbank
    records

    Attributes:
        parent_gbk: GBK
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
    """

    # TODO: replace any with GBK after restructuring
    def __init__(self):
        self.parent_gbk: Optional[Any] = None
        # contig edge is optional, proto_core does not have it
        self.contig_edge: Optional[bool] = None
        self.nt_start: Optional[int] = None
        self.nt_stop: Optional[int] = None
        self.product: Optional[str] = None

    def save(self, type: str, commit=True):
        """Stores this BGCRecord in the database

        Args:
            type (str):  type of antiSMASH record
            commit (bool, optional): commit immediately after executing the insert
            query. Defaults to True.
        """

        bgc_record_table = DB.metadata.tables["bgc_record"]

        contig_edge = None
        if self.contig_edge is not None:
            contig_edge = self.contig_edge

        parent_gbk_id = None
        if self.parent_gbk is not None and self.parent_gbk._db_id is not None:
            parent_gbk_id = self.parent_gbk._db_id

        insert_query = (
            bgc_record_table.insert()
            .prefix_with("OR REPLACE")
            .values(
                gbk_id=parent_gbk_id,
                contig_edge=contig_edge,
                nt_start=self.nt_start,
                nt_stop=self.nt_stop,
                type=type,
            )
            .compile()
        )

        DB.execute(insert_query)

        if commit:
            DB.commit()

    def parse_bgc_record(self, feature: SeqFeature, parent_gbk: Optional[Any]):
        """Parses a BGC record locale info

        Args:
            feature (SeqFeature): SeqFeature, any type BGC record

        Raises:
            InvalidGBKError: Invalid or missing fields in SeqFeature
        """

        if "contig_edge" in feature.qualifiers:
            contig_edge_qualifier = feature.qualifiers["contig_edge"][0]

            if contig_edge_qualifier == "True":
                self.contig_edge = True
            else:
                self.contig_edge = False

        self.nt_start = feature.location.start
        self.nt_stop = feature.location.end

        if "product" not in feature.qualifiers:
            logging.error("product qualifier not found in feature!")
            raise InvalidGBKError()

        self.product = feature.qualifiers["product"][0]

        # add parent gbk if available
        if parent_gbk is not None:
            self.parent_gbk = parent_gbk
