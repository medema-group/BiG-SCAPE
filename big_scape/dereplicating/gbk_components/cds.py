"""Contains code relevant to the CDS component"""

# from python

# from dependencies
from Bio.Seq import Seq

# from other modules

# from this module


class CDS:
    """
    Class to describe a CDS component

    Attributes:
        nt_start (int): start position of the CDS in the nucleotide sequence
        nt_end (int): end position of the CDS in the nucleotide sequence
        biosynthetic (bool): whether the CDS is core biosynthetic
        aa_sequence (Bio.Seq): amino acid sequence of the CDS

    """

    def __init__(
        self, nt_start: int, nt_end: int, biosynthetic: bool, aa_sequence: Seq
    ) -> None:
        self.nt_start: int = nt_start
        self.nt_end: int = nt_end
        self.biosynthetic: bool = biosynthetic
        self.aa_sequence: Seq = aa_sequence
