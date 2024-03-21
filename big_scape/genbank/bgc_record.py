"""Contains a class for a BGC record, a common set of qualifiers/attributes amongst the
AntiSMASh genbank records
"""

# from python
from __future__ import annotations
from typing import Optional, Sequence, TYPE_CHECKING
import logging

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from big_scape.cli.constants import ANTISMASH_CLASSES
from big_scape.data import DB
from big_scape.errors import InvalidGBKError
from big_scape.genbank.cds import CDS

import big_scape.enums as bs_enums

# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import (
        GBK,
        Region,
        ProtoCluster,
        ProtoCore,
    )  # imported earlier in file_input.load_files
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
        merged: bool
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
        self.merged: bool = False

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

    def get_hsps(self, return_all=False) -> list[HSP]:
        """Get a list of all hsps in this region

        Args:
            return_all (bool): If set to true, returns all HSP regardless of coordinate
            information. Defaults to False

        Returns:
            list[HSP]: List of all hsps in this region
        """
        domains: list[HSP] = []
        for cds in self.get_cds_with_domains(return_all=return_all):
            if len(cds.hsps) > 0:
                if cds.strand == 1:
                    domains.extend(cds.hsps)
                elif cds.strand == -1:
                    domains.extend(cds.hsps[::-1])
        return domains

    def get_cds_start_stop(self) -> tuple[int, int]:
        """Get cds ORF number of record start and stop with respect to full region

        Obtained cds slice starts counting at one and is __inclusive__

        Args:
            record (BGCRecord): record to find bounds for

        Returns:
            tuple[int, int]: start and stop of record in cds number
        """
        if self.parent_gbk is None:
            raise AttributeError("Record parent GBK is not set!")

        gbk = self.parent_gbk
        all_cds = gbk.genes

        record_start = 1
        record_stop = len(all_cds)
        # check if record contains all cds
        if all_cds[0].nt_start >= self.nt_start and all_cds[-1].nt_stop <= self.nt_stop:
            return record_start, record_stop

        for idx, cds in enumerate(all_cds):
            if cds.nt_start < self.nt_start:
                record_start = idx + 2
            if cds.nt_stop > self.nt_stop:
                record_stop = idx
                break
        return record_start, record_stop

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

        if not hasattr(self, "category"):
            self.category: Optional[str] = None

        # why
        if hasattr(self, "merged_number"):
            number = self.merged_number
        else:
            number = self.number

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
                record_number=number,
                contig_edge=contig_edge,
                nt_start=self.nt_start,
                nt_stop=self.nt_stop,
                product=self.product,
                category=self.category,
                record_type=record_type,
                merged=self.merged,
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
    def parse_products(feature: SeqFeature) -> str:
        """Parse a set of products from a BGC record from antismash V5 and above.

        Args:
            feature (SeqFeature): BGC record feature to parse

        Returns:
            str: Singular string representing the product type
        """

        if "product" not in feature.qualifiers:
            logging.error("product qualifier not found in feature!")
            raise InvalidGBKError()

        products: list[str] = feature.qualifiers["product"]

        # single product? just return it
        if len(products) == 1:
            return products.pop()

        products.sort()

        # TODO: check the status here
        # return easy hybrids
        # if "other" not in products:
        return ".".join(products)

        # TODO: clean up
        # in all other cases we have an 'other' classification. for the rest of the
        # cases we can remove that and just parse the products again
        # products.remove("other")
        # return BGCRecord.parse_products(products)

    @staticmethod
    def parse_products_as4(feature: SeqFeature) -> str:
        """Parse a set of products from a BGC record from antismash V4 and under.

        Args:
            feature (SeqFeature): BGC record feature to parse

        Returns:
            str: Singular string representing the product type
        """

        if "product" not in feature.qualifiers:
            logging.error("product qualifier not found in feature!")
            raise InvalidGBKError()

        as4_products = []
        for product_type in ANTISMASH_CLASSES:
            for product in ANTISMASH_CLASSES[product_type]:
                as4_products.append(product)

        product = feature.qualifiers["product"][0]

        if "-" not in product:
            return product

        elif product in as4_products:
            return product

        else:
            products = []
            for as4_product in as4_products:
                if as4_product in product:
                    products.append(as4_product)
            products.sort()
            return ".".join(products)

    @staticmethod
    def parse_common(
        feature: SeqFeature,
    ) -> tuple[int, int, Optional[bool]]:
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

        return nt_start, nt_stop, contig_edge


def get_sub_records(
    region: Region, record_type: bs_enums.genbank.RECORD_TYPE
) -> Sequence[BGCRecord]:
    """Return a list of BGCRecords of the specified type from the given region

    This function recurses through a record, returning any record of a given type. If
    the record is of the specified type, it is returned. If it is not, the function
    recurses through the record's children, returning any records of the specified type.

    If a given type of record is not found, the function will return the region instead

    Args:
        region (Region): region to search
        type (bs_enums.genbank.RECORD_TYPE): type of record to return

    Returns:
        list[BGCRecord]: list of records of the specified type
    """
    if region is None:
        return []

    if record_type == bs_enums.genbank.RECORD_TYPE.REGION:
        return [region]

    if region.cand_clusters is None:
        return [region]

    cand_clusters = [
        cand_cluster
        for cand_cluster in region.cand_clusters.values()
        if cand_cluster is not None
    ]

    if len(cand_clusters) == 0:
        return [region]

    # TODO: check if properly implemented, or at all
    #  use same strategy as for region if kind is neighbouring (several cores in a record)
    # and protocluster if kind is single, interleaved or chemical hybrid (one core per record)
    if record_type == bs_enums.genbank.RECORD_TYPE.CAND_CLUSTER:
        return cand_clusters

    # for protoclusters and protocores we need to keep track of the numbers, since they may
    # be repeated in different candidate clusters

    proto_clusters: list[ProtoCluster] = []
    proto_cluster_numbers: list[int] = []
    for cand_cluster in cand_clusters:
        if cand_cluster.proto_clusters is None:
            continue

        for proto_cluster in cand_cluster.proto_clusters.values():
            if (
                proto_cluster is not None
                and proto_cluster.number not in proto_cluster_numbers
            ):
                proto_clusters.append(proto_cluster)
                proto_cluster_numbers.append(proto_cluster.number)

    if len(proto_clusters) == 0:
        return cand_clusters

    if record_type == bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER:
        return proto_clusters

    proto_cores: list[ProtoCore] = []
    proto_core_numbers: list[int] = []
    for proto_cluster in proto_clusters:
        if proto_cluster.proto_core is None:
            continue
        for proto_core in proto_cluster.proto_core.values():
            if proto_core is not None and proto_core.number not in proto_core_numbers:
                proto_cores.append(proto_core)
                proto_core_numbers.append(proto_core.number)

    if len(proto_cores) == 0:
        return proto_clusters

    return proto_cores
