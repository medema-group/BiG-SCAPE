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
from big_scape.data import DB
from big_scape.errors import InvalidGBKError
from big_scape.genbank.cds import CDS

# from this module


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files
    from big_scape.hmm import HSP  # imported earlier in genbank.CDS


class BGCRecord:
    """Describes a common set of qualifiers/attributes among all ඞ AntiSMASH genbank
    records

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
        self.parent_gbk = parent_gbk
        self.number = number
        # contig edge is optional, proto_core does not have it
        self.contig_edge = contig_edge
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.product = product

        # for database operations
        self._db_id: Optional[int] = None

        # for networking
        self._families: dict[float, int] = {}

    def get_cds(self, return_all=False, reverse=False) -> list[CDS]:
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
            nt_start and nt_stop or all CDS if return_all is true
        """
        if self.parent_gbk is None:
            raise ValueError("BGCRegion does not have a parent")

        parent_gbk_cds: list[CDS] = self.parent_gbk.genes

        if return_all:
            # TODO: I don't like this solution. maybe go back to the more difficult one
            if reverse:
                return list(reverse(self.parent_gbk.genes))

            return list(self.parent_gbk.genes)

        if self.nt_start is None or self.nt_stop is None:
            raise ValueError("Cannot CDS from region with no position information")

        record_cds: list[CDS] = []

        if reverse:
            step = -1
        else:
            step = 1

        for cds in parent_gbk_cds[::step]:
            if cds.nt_start < self.nt_start:
                continue

            if cds.nt_stop > self.nt_stop:
                continue

            record_cds.append(cds)

        return record_cds

    def get_cds_with_domains(self, return_all=False, reverse=False) -> list[CDS]:
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
            nt_start and nt_stop or all CDS if return_all is true
        """
        if self.parent_gbk is None:
            raise ValueError("BGCRegion does not have a parent")

        parent_gbk_cds: list[CDS] = self.parent_gbk.genes

        if return_all:
            # TODO: I don't like this solution. maybe go back to the more difficult one
            if reverse:
                return [
                    cds for cds in reversed(self.parent_gbk.genes) if len(cds.hsps) > 0
                ]

            return [cds for cds in self.parent_gbk.genes if len(cds.hsps) > 0]

        if self.nt_start is None or self.nt_stop is None:
            raise ValueError("Cannot CDS from region with no position information")

        record_cds: list[CDS] = []

        if reverse:
            step = -1
        else:
            step = 1

        for cds in parent_gbk_cds[::step]:
            if len(cds.hsps) == 0:
                continue

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
        domains: list[HSP] = []
        for cds in self.get_cds_with_domains():
            if len(cds.hsps) > 0:
                domains.extend(cds.hsps)
        return domains

    def save_record(
        self, record_type: str, parent_id: Optional[int] = None, commit=True
    ) -> None:
        """Stores this BGCRecord in the database

        Args:
            type (str):  type of antiSMASH record
            parent_id (Optional[int]): database ID of the parent bgc record. None if
            this record is a region. Defaults to None.
            commit (bool, optional): commit immediately after executing the insert
            query. Defaults to True.
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        bgc_record_table = DB.metadata.tables["bgc_record"]

        contig_edge = None
        if self.contig_edge is not None:
            contig_edge = self.contig_edge

        gbk_id = None
        if self.parent_gbk is not None and self.parent_gbk._db_id is not None:
            gbk_id = self.parent_gbk._db_id

        insert_query = (
            bgc_record_table.insert()
            .prefix_with("OR REPLACE")
            .values(
                gbk_id=gbk_id,
                parent_id=parent_id,
                record_number=self.number,
                contig_edge=contig_edge,
                nt_start=self.nt_start,
                nt_stop=self.nt_stop,
                product=self.product,
                record_type=record_type,
            )
            .returning(bgc_record_table.c.id)
            .compile()
        )

        cursor_result = DB.execute(insert_query, False)

        # get return value
        return_row = cursor_result.fetchone()

        if return_row is None:
            raise RuntimeError("No return value from insert query")

        self._db_id = return_row[0]

        if commit:
            DB.commit()

    def get_attr_dict(self) -> dict[str, object]:
        """Gets a dictionary of attributes, used for adding to network nodes later

        Returns:
            dict[str, object]: dictionary of attributes
        """
        attr_dict: dict[str, object] = {
            "product": self.product,
            "contig_edge": self.contig_edge,
        }

        for cutoff, label in self._families.items():
            attr_dict[str(cutoff)] = label

        return attr_dict

    def __repr__(self) -> str:
        return f"{self.parent_gbk} Record (superclass) {self.nt_start}-{self.nt_stop}"

    def __hash__(self, record_type="BGCRecord") -> int:
        return hash(
            (
                self.parent_gbk,
                record_type,
                self.nt_start,
                self.nt_stop,
            )
        )

    @staticmethod
    def parse_products(products: set[str]) -> str:
        """Parse a set of products from a BGC record. Used in cases where there are hybrid products

        Args:
            products (set[str]): set of products

        Returns:
            str: Singular string representing the product type
        """
        # single product? just return it
        if len(products) == 1:
            return products.pop()

        # return easy hybrids
        if "other" not in products:
            return ".".join(products)

        # in all other cases we have an 'other' classification. for the rest of the
        # cases we can remove that and just parse the products again
        products.remove("other")
        return BGCRecord.parse_products(products)

    @staticmethod
    def parse_common(
        feature: SeqFeature,
    ) -> tuple[int, int, Optional[bool], str]:
        """Parse and return the common attributes of BGC records in a GBK feature

        Args:
            feature (SeqFeature): Seq feature of a BGC record to parse

        Raises:
            InvalidGBKError: Raised when a field that was expected to be present in the
            feature is not present

        Returns:
            tuple[int, int, Optional[bool], str]: nt_start, nt_stop, contig_edge,
            product
        """
        nt_start = feature.location.start
        nt_stop = feature.location.end

        contig_edge = None
        if "contig_edge" in feature.qualifiers:
            contig_edge_qualifier = feature.qualifiers["contig_edge"][0]

            contig_edge = contig_edge_qualifier == "True"

        if "product" not in feature.qualifiers:
            logging.error("product qualifier not found in feature!")
            raise InvalidGBKError()

        # record may have multiple products. handle them here
        products = set(feature.qualifiers["product"][0].split("-"))

        product = BGCRecord.parse_products(products)

        return nt_start, nt_stop, contig_edge, product