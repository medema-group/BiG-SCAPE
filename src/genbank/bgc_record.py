"""Contains a class for a BGC record, a common set of qualifiers/attributes amongst the
AntiSMASh genbank records
"""

# from python
from __future__ import annotations
from typing import Optional, TYPE_CHECKING
import logging

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB
from src.errors import InvalidGBKError
from src.genbank.cds import CDS

# from this module


# from circular imports
if TYPE_CHECKING:
    from src.genbank import GBK  # imported earlier in src.file_input.load_files
    from src.hmm import HSP  # imported earlier in src.genbank.CDS


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

    def __init__(self):
        self.parent_gbk: Optional[GBK] = None
        # contig edge is optional, proto_core does not have it
        self.contig_edge: Optional[bool] = None
        self.nt_start: Optional[int] = None
        self.nt_stop: Optional[int] = None
        self.product: Optional[str] = None

    def get_cds(self, return_all=False) -> list[CDS]:
        """Get a list of CDS that lie within the coordinates specified in this region
        from the parent GBK class

        Args:
            return_all (bool): If set to true, returns all CDS regardless of coordinate
            information. Defaults to False

        Raises:
            ValueError: Raised if this class contains no position information or if this
            record does not have a parent

        Returns:
            list[CDS]: A list of CDS that lie only within the coordinates specified by
            nt_start and nt_stop
        """
        if self.parent_gbk is None:
            raise ValueError("BGCRegion does not have a parent")

        parent_gbk_cds: list[CDS] = self.parent_gbk.genes

        if return_all:
            return self.parent_gbk.genes

        if self.nt_start is None or self.nt_stop is None:
            raise ValueError("Cannot CDS from region with no position information")

        record_cds = []
        for cds in parent_gbk_cds:
            if cds.nt_start < self.nt_start:
                continue

            if cds.nt_stop > self.nt_stop:
                continue

            record_cds.append(cds)

        return record_cds

    def get_hsps(self) -> list[HSP]:
        """Get a list of all hsps in this region

        Returns:
            list[HSP]: List of all hsps in this region
        """
        domains = []
        for cds in self.get_cds():
            domains.extend(cds.hsps)
        return domains

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

    def parse_bgc_record(self, feature: SeqFeature, parent_gbk: Optional[GBK]):
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

    def __repr__(self) -> str:
        return f"{self.parent_gbk} Record (superclass) {self.nt_start}-{self.nt_stop} "

    def __hash__(self, record_type="BGCRecord") -> int:
        # return a hash of a tuple containing identifying properties
        return hash(
            (
                self.parent_gbk,
                record_type,
                self.nt_start,
                self.nt_stop,
            )
        )
