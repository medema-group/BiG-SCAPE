""" Module containing code to load and store antiSMASH CDSs"""

# from python
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
        domains: list[str]
    """

    def __init__(self, nt_start: int, nt_stop: int):
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.gene_kind: Optional[str] = None
        self.strand: Optional[int] = None
        self.aa_seq: Optional[str] = None
        self.domains: list = []

    # def add_domain(self, domain: Domain):
    #     """Add a domain object to this CDS"""

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
        nt_end = int(feature.location.end)
        strand = int(feature.strand)

        cds = cls(nt_start, nt_end)
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
