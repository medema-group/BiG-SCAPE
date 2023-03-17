"""Contains a class for a BGC record, a common set of qualifiers/attributes amongst the AntiSMASH
genbank records
"""

# from python
from __future__ import annotations
from typing import Optional

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.data import DB

# from this module


class BGCRecord:
    """Describes a common set of qualifiers/attributes amongst all AntiSMASH genbank records"""

    def __init__(self):
        self.contig_edge: Optional[bool] = None
        self.nt_start: Optional[int] = None
        self.nt_stop: Optional[int] = None

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

    def parse_location(self, feature: SeqFeature):
        """Parses a BGC record locale info"""

        if "contig_edge" in feature.qualifiers:
            contig_edge_qualifier = feature.qualifiers["contig_edge"][0]

            if contig_edge_qualifier == "True":
                self.contig_edge = True
            else:
                self.contig_edge = False

        self.nt_start = feature.location.start
        self.nt_stop = feature.location.end
