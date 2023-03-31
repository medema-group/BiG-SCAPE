""" Module containing code to load and store antiSMASH CDSs"""

# from python
from __future__ import annotations
from itertools import combinations
import logging
from typing import Optional

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError


class CDS:
    """
    Class to describe a CDS within an antiSMASH GBK

    Attributes:
        gene_kind: str
        strand: Bool
        nt_start: int
        nt_sop: int
        aa_seq: SeqRecord.seq
        hsps: list[HSP]
    """

    def __init__(self, nt_start: int, nt_stop: int):
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.gene_kind: Optional[str] = None
        self.strand: Optional[int] = None
        self.aa_seq: str = ""
        self.hsps: list = []

    @classmethod
    def parse(cls, feature: SeqFeature):
        """Creates a cds object from a region feature in a GBK file"""

        if feature.type != "CDS":
            logging.error(
                "Feature is not of correct type! (expected: region, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        nt_start = int(feature.location.start)
        nt_stop = int(feature.location.end)
        strand = int(feature.location.strand)

        cds = cls(nt_start, nt_stop)
        cds.strand = strand

        if "translation" not in feature.qualifiers:
            logging.error("translation qualifier not found in cds feature!")
            raise InvalidGBKError()

        aa_seq = str(feature.qualifiers["translation"][0])
        cds.aa_seq = aa_seq

        if "gene_kind" in feature.qualifiers:
            gene_kind = str(feature.qualifiers["gene_kind"][0])
            cds.gene_kind = gene_kind

        return cds

    @staticmethod
    def has_overlap(cds_a: CDS, cds_b: CDS) -> bool:
        """Return True if there is overlap between this

        Args:
            cds_b (CDS): CDS to compare

        Returns:
            bool: whether there is overlap between this cds and another
        """
        has_overlap = CDS.len_overlap(cds_a, cds_b) > 0
        return has_overlap

    @staticmethod
    def len_overlap(cds_a: CDS, cds_b: CDS) -> int:
        """Return the length of the overlap between this CDS and another

        Args:
            cds_b (CDS): CDS to compare

        Returns:
            int: length of the overlap between this CDS and another
        """

        if cds_a.nt_start < cds_b.nt_start:
            left = cds_b.nt_start
        else:
            left = cds_a.nt_start

        if cds_a.nt_stop > cds_b.nt_stop:
            right = cds_b.nt_stop
        else:
            right = cds_a.nt_stop

        # limit to > 0
        return max(0, right - left)

    @staticmethod
    def filter_overlap(cds_list: list[CDS]):
        # TODO: document

        # working with lists here is kind of iffy. in this case we are keeping track of\
        # which CDS we want to remove from the original list later on
        del_list = set()
        # find all combinations of cds to check for overlap
        cds_a: CDS
        cds_b: CDS
        for cds_a, cds_b in combinations(cds_list, 2):
            a_len = len(cds_a.aa_seq)
            b_len = len(cds_b.aa_seq)
            shortest_len = min(a_len, b_len)

            # do not add to remove list if there is no overlap at all
            if not CDS.has_overlap(cds_a, cds_b):
                continue

            # calculate overlap
            nt_overlap = CDS.len_overlap(cds_a, cds_b)
            aa_overlap = nt_overlap / 3

            # allow the overlap to be as large as 10% of the shortest CDS.
            if aa_overlap > 0.1 * shortest_len:
                if a_len > b_len:
                    del_list.add(cds_b)
                else:
                    del_list.add(cds_a)

        # remove any entries that need to be removed
        for cds in del_list:
            cds_list.remove(cds)

        return cds_list
