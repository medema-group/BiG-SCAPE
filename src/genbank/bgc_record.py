"""Contains a class for a BGC record, a common set of qualifiers/attributes amongst the AntiSMASH
genbank records
"""

# from python
from __future__ import annotations
import logging

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB
from src.errors import InvalidGBKError

# from this module


class BGCRecord:
    """Describes a common set of qualifiers/attributes amongst all AntiSMASH genbank records"""

    def __init__(self, contig_edge: bool, nt_start: int, nt_stop: int, product: str):
        self.contig_edge = contig_edge
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.product = product

    def save(self, type: str, commit=True):
        """Stores this BGCRecord in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        bgc_record_table = DB.metadata.tables["bgc_record"]

        contig_edge = None
        if self.contig_edge is not None:
            contig_edge = self.contig_edge

        insert_query = (
            bgc_record_table.insert()
            .prefix_with("OR REPLACE")
            .values(
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

    @staticmethod
    def parse_bgc_record(feature: SeqFeature):
        """Parses a BGC record locale info"""

        if "contig_edge" in feature.qualifiers:
            contig_edge_qualifier = feature.qualifiers["contig_edge"][0]

            if contig_edge_qualifier == "True":
                contig_edge = True
            else:
                contig_edge = False

        nt_start = feature.location.start
        nt_stop = feature.location.end

        if "product" not in feature.qualifiers:
            logging.error("product qualifier not found in feature!")
            raise InvalidGBKError()

        product = feature.qualifiers["product"][0]

        return (contig_edge, nt_start, nt_stop, product)
