""" Module containing code to load and store antiSMASH CDSs"""

# from python
from __future__ import annotations
from itertools import combinations
import logging
from typing import Optional, TYPE_CHECKING

# from dependencies
from Bio.Data.CodonTable import TranslationError
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq

# from other modules
from src.errors import InvalidGBKError
from src.data import DB


# from circular imports
if TYPE_CHECKING:
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


class CDS:
    """
    Class to describe a CDS within an antiSMASH GBK

    Attributes:
        gene_kind: str
        strand: Bool
        nt_start: int
        nt_stop: int
        aa_seq: SeqRecord.seq
        hsps: list[HSP]
    """

    def __init__(self, nt_start: int, nt_stop: int):
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.parent_gbk: Optional[GBK] = None
        self.gene_kind: Optional[str] = None
        self.strand: Optional[int] = None
        self.aa_seq: str = ""
        self.hsps: list = []

        # db specific fields
        self._db_id: Optional[int] = None

    def save(self, commit=True):
        """Saves this CDS to the database and optionally executes a commit

        Args:
            commit (bool, optional): Whether to commit immediately after inserting this
            CDS. Defaults to True.
        """

        # get parent gbk id if available
        parent_gbk_id = None
        if self.parent_gbk is not None and self.parent_gbk._db_id is not None:
            parent_gbk_id = self.parent_gbk._db_id

        cds_table = DB.metadata.tables["cds"]
        insert_query = (
            cds_table.insert()
            .returning(cds_table.c.id)
            .values(
                gbk_id=parent_gbk_id,
                nt_start=self.nt_start,
                nt_stop=self.nt_stop,
                strand=self.strand,
                gene_kind=self.gene_kind,
                aa_seq=self.aa_seq,
            )
            .compile()
        )

        # in the above query we add a returning statement. This makes it so that the
        # sqlite engine will be in the middle of a transaction (trying to give us the
        # returned value) and means we cannot commit. We will do this after digesting
        # the reutrn statement further on
        cursor_result = DB.execute(insert_query, False)

        # get return value
        return_row = cursor_result.fetchone()
        self._db_id = return_row[0]

        # only now that we have handled the return we can commit
        if commit:
            DB.commit()

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: GBK):
        """Creates a cds object from a region feature in a GBK file"""

        if feature.type != "CDS":
            logging.error(
                "Feature is not of correct type! (expected: region, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        # if fuzzy positions present follows behaviour described
        # here: https://biopython.org/docs/1.81/api/Bio.SeqFeature.html
        nt_start = int(feature.location.start)
        nt_stop = int(feature.location.end)
        strand = int(feature.location.strand)

        cds = cls(nt_start, nt_stop)
        cds.strand = strand

        nt_seq = feature.location.extract(parent_gbk.nt_seq)

        # if translation exists -> check if there is a match with biopython translation
        # if no match, pick AS version and issue warning

        if "translation" in feature.qualifiers:
            aa_seq = str(feature.qualifiers["translation"][0])

            # translate own
            transl_nt_seq = translate(feature, nt_seq)

            # TODO make better warnings
            if not transl_nt_seq:
                logging.warning("could not translate CDS nucleotide sequence")

            # check if both translations match
            if transl_nt_seq and not str(nt_seq) == str(transl_nt_seq):
                logging.warning(
                    "translation of CDS sequence from antiSMASH"
                    " and biopython do not match"
                )

            cds.aa_seq = aa_seq

        else:
            # if nt_seq is not divisible by 3, we discard this CDS
            transl_nt_seq = translate(feature, nt_seq)
            if not transl_nt_seq:
                logging.warning(
                    "translation qualifier not found in cds feature"
                    " and could not be computed, therefore this CDS"
                    " feature is being discarded"
                )
                return None

        if "gene_kind" in feature.qualifiers:
            gene_kind = str(feature.qualifiers["gene_kind"][0])
            cds.gene_kind = gene_kind

        # add parent if it exists
        if parent_gbk is not None:
            cds.parent_gbk = parent_gbk

        return cds

    @staticmethod
    def has_overlap(cds_a: CDS, cds_b: CDS) -> bool:
        """Return True if there is overlap between this

        Args:
            cds_b (CDS): CDS to compare

        Returns:
            bool: whether there is overlap between this cds and another
        """
        has_overlap = CDS.len_nt_overlap(cds_a, cds_b) > 0
        return has_overlap

    @staticmethod
    def len_nt_overlap(cds_a: CDS, cds_b: CDS) -> int:
        """Return the length of the nucleotide seq overlap between two CDSs

        Args:
            cds_a,b (CDS): CDSs to compare

        Returns:
            int: length of the overlap between this CDS and another
        """
        if cds_a.strand != cds_b.strand:
            return 0

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
    def filter_overlap(cds_list: list[CDS], perc_overlap: float):
        """From a list of CDSs, filter out overlapping CDSs if overlap len exceeds
        a percentage of the shortest total CDS len

        Args:
            cds_list (list[CDS]): list of CDS
            perc_overlap (float): threshold to filter

        Returns:
            _type_: list[CDS]
        """

        # working with lists here is kind of iffy. in this case we are keeping track of\
        # which CDS we want to remove from the original list later on
        del_list = set()
        # find all combinations of cds to check for overlap
        cds_a: CDS
        cds_b: CDS
        for cds_a, cds_b in combinations(cds_list, 2):
            a_aa_len = len(cds_a.aa_seq)
            b_aa_len = len(cds_b.aa_seq)
            shortest_aa_len = min(a_aa_len, b_aa_len)

            # do not add to remove list if cds are on a different strand
            if cds_a.strand != cds_b.strand:
                continue

            # do not add to remove list if there is no overlap at all
            if not CDS.has_overlap(cds_a, cds_b):
                continue

            # calculate overlap in nt and convert to aa
            nt_overlap = CDS.len_nt_overlap(cds_a, cds_b)
            aa_overlap = nt_overlap / 3

            # allow the aa overlap to be as large as X% of the shortest CDS.
            if aa_overlap > perc_overlap * shortest_aa_len:
                if a_aa_len > b_aa_len:
                    del_list.add(cds_b)
                else:
                    del_list.add(cds_a)

        # remove any entries that need to be removed
        for cds in del_list:
            cds_list.remove(cds)

        return cds_list


def translate(feature: SeqFeature, nt_seq: Seq):
    """Translate a CDS sequence if CDS is valid

    cds = True: checks the sequence starts with a valid alternative start codon
    (which will be translated as methionine, M), that the sequence length is a
    multiple of three, and that there is a single in frame stop codon at the end
    (this will be excluded from the protein sequence, regardless of the to_stop option).
    If these tests fail, an exception is raised.
    (https://biopython.org/docs/1.81/api/Bio.Seq.html)


    Args:
        feature (SeqFeature): CDS SeqFeature
        nt_seq (Seq): CDS nucleotide sequence

    Returns:
        _type_: Seq or None
    """

    transl_table = None

    if "transl_table" in feature.qualifiers.keys():
        transl_table = feature.qualifiers["transl_table"][0]

    try:
        transl_nt_seq = nt_seq.translate(table=transl_table, cds=True)
    except TranslationError:
        transl_nt_seq = None

    return transl_nt_seq
