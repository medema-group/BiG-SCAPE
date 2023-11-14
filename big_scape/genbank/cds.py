""" Module containing code to load and store antiSMASH CDSs"""

# from python
from __future__ import annotations
import logging
from typing import Optional, TYPE_CHECKING
import warnings

# from dependencies
from Bio.SeqFeature import SeqFeature
from Bio.Seq import Seq
from Bio import BiopythonWarning

# from other modules
from big_scape.errors import InvalidGBKError
from big_scape.data import DB
from big_scape.hmm import HSP
from big_scape.enums import SOURCE_TYPE

# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files


class CDS:
    """
    Class to describe a CDS within an antiSMASH GBK

    Attributes:
        nt_start: int
        nt_stop: int
        orf_num: int
        parent_gbk: GBK
        gene_kind: str
        strand: Bool
        aa_seq: SeqRecord.seq
        hsps: list[HSP]
    """

    def __init__(self, nt_start: int, nt_stop: int):
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.orf_num: Optional[int] = None
        self.parent_gbk: Optional[GBK] = None
        self.gene_kind: Optional[str] = None
        self.strand: Optional[int] = None
        self.aa_seq: str = ""
        self.hsps: list[HSP] = []

        # db specific fields
        self._db_id: Optional[int] = None

    def add_hsp_overlap_filter(self, new_hsp: HSP, domain_overlap_cutoff=0.1) -> None:
        """Adds a HSP to this CDS. Performs overlap cutoff filtering by calculating the
        percentage overlap of the incoming HSP with other HSPs in this CDS.

        If the percentage overlap is greater than the cutoff, this keeps whichever HSP
        has the higher score. If scores are equal, this keeps the HSP with the earliest
        start position. Bitscores are rounded to 1 decimal position when compared.

        The above behavior should mirror BiG-SCAPE 1.0 behavior

        Args:
            hsp (hmmer.hsp): The HSP to be added to this CDS
            overlap_cutoff (float, optional): cutoff threshold for overlap. Defaults to
            0.1

        """
        # if no hsps added yet, just add and continue
        if len(self.hsps) == 0:
            self.hsps.append(new_hsp)
            return

        delete_list = []

        for hsp_idx, old_hsp in enumerate(self.hsps):
            # just add it if there is no overlap at all
            if not HSP.has_overlap(old_hsp, new_hsp):
                continue

            # there is some overlap. calculate how much
            overlap_aa = HSP.len_overlap(old_hsp, new_hsp)

            overlap_perc_old = overlap_aa / (old_hsp.env_stop - old_hsp.env_start)
            overlap_perc_new = overlap_aa / (new_hsp.env_stop - new_hsp.env_start)

            # neither over cutoff?
            if (
                overlap_perc_old < domain_overlap_cutoff
                and overlap_perc_new < domain_overlap_cutoff
            ):
                continue

            score_old = round(old_hsp.score, 1)
            score_new = round(new_hsp.score, 1)

            # keep old if it has a better score
            if score_old > score_new:
                return

            # replace old if new is better
            if score_new > score_old:
                delete_list.append(hsp_idx)
                continue

            # if scores are equal, keep the one with the lower nt_start
            if new_hsp.env_start < old_hsp.env_start:
                delete_list.append(hsp_idx)
                continue

        # go through this in reverse order otherwise we mess everything up
        for deleted_idx in delete_list[::-1]:
            del self.hsps[deleted_idx]

        # if we got through all of that without the function, we never replaced an HSP
        # so add a new one here
        self.hsps.append(new_hsp)

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
                orf_num=self.orf_num,
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

        if return_row is None:
            raise RuntimeError("No return value from insert query")

        self._db_id = return_row[0]

        # only now that we have handled the return we can commit
        if commit:
            DB.commit()

    def __eq__(self, __o) -> bool:
        if not isinstance(__o, CDS):
            raise NotImplementedError()

        # exception: no hsps in both
        # TODO: time to stop abusing dunder methods. use custom sorted lists
        if len(self.hsps) + len(__o.hsps) == 0:
            return self.parent_gbk == __o.parent_gbk and self.orf_num == __o.orf_num

        return self.hsps == __o.hsps

    def __gt__(self, __o) -> bool:
        if not isinstance(__o, CDS):
            raise NotImplementedError()

        return self.nt_start > __o.nt_start

    def __hash__(self) -> int:
        return hash(tuple(self.hsps))

    def __repr__(self) -> str:
        if self.parent_gbk is None:
            parent_gbk_str = "ORPHAN"
        else:
            parent_gbk_str = str(self.parent_gbk.path.name)

        domain_list = " ".join([f"{hsp.domain[2:]: <8}" for hsp in self.hsps])

        return (
            f"{parent_gbk_str}_CDS{self.orf_num}, {self.nt_start}-{self.nt_stop}:"
            f"{'+' if self.strand == 1 else '-'} {self.gene_kind: <12} - {domain_list}"
        )

    @classmethod
    def parse(
        cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None
    ) -> Optional[CDS]:
        """Creates a cds object from a region feature in a GBK file

        Note that this will not add biosynthetic information if no parent gbk is passed
        """

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

        # add parent if it exists
        if parent_gbk is None:
            return cds

        cds.parent_gbk = parent_gbk
        cds.gene_kind = ""
        if parent_gbk.as_version == "4":
            if "sec_met" in feature.qualifiers:
                for sec_met_value in feature.qualifiers["sec_met"]:
                    if "Kind" in sec_met_value:
                        cds.gene_kind = sec_met_value[6:]  # trim "Kind: "
                        break

        if "gene_kind" in feature.qualifiers:
            cds.gene_kind = str(feature.qualifiers["gene_kind"][0])

        nt_seq = feature.location.extract(parent_gbk.nt_seq)

        # if translation exists -> check if there is a match with biopython translation
        # if no match, pick AS version and issue warning

        if "translation" in feature.qualifiers:
            aa_seq = str(feature.qualifiers["translation"][0])

            correct_translation = check_translation(aa_seq, nt_seq, feature)

            if not correct_translation:
                correct_skip_start_translation = check_translation(
                    aa_seq[1:], nt_seq[1:], feature
                )

                if (
                    not correct_skip_start_translation
                    and not parent_gbk.source_type == SOURCE_TYPE.MIBIG
                ):
                    logging.warning(
                        "CDS (%s, %s) from %s:"
                        " translation provided by antiSMASH and generated by biopython"
                        " do not match, consider checking if there is something"
                        " special with this CDS",
                        nt_start,
                        nt_stop,
                        parent_gbk.path,
                    )

            cds.aa_seq = aa_seq
            return cds

        # if not translation available -> we try to generate it

        transl_nt_seq = get_translation(feature, nt_seq)

        # we can't work with a CDS that has no AA seq so ignore it by returning None
        if transl_nt_seq is None:
            logging.warning(
                "CDS (%s, %s) from %s:"
                " translation not found in cds feature and could not be generated,"
                " therefore this CDS feature is being discarded",
                nt_start,
                nt_stop,
                parent_gbk.path,
            )
            return None

        cds.aa_seq = transl_nt_seq
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
    def load_all(gbk_dict: dict[int, GBK]) -> None:
        """Load all Region objects from the database

        This function populates the region objects in the GBKs provided in the input
        gbk_dict

        Args:
            region_dict (dict[int, GBK]): Dictionary of Region objects with database ids
            as keys. Used for parenting
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        cds_table = DB.metadata.tables["cds"]

        region_select_query = (
            cds_table.select()
            .add_columns(
                cds_table.c.id,
                cds_table.c.gbk_id,
                cds_table.c.nt_start,
                cds_table.c.nt_stop,
                cds_table.c.orf_num,
                cds_table.c.strand,
                cds_table.c.gene_kind,
                cds_table.c.aa_seq,
            )
            .order_by(cds_table.c.orf_num)
            .where(cds_table.c.gbk_id.in_(gbk_dict.keys()))
            .compile()
        )

        cursor_result = DB.execute(region_select_query)

        for result in cursor_result.all():
            new_cds = CDS(result.nt_start, result.nt_stop)

            new_cds._db_id = result.id

            new_cds.parent_gbk = gbk_dict[result.gbk_id]
            new_cds.strand = result.strand
            new_cds.aa_seq = result.aa_seq
            new_cds.orf_num = result.orf_num
            new_cds.gene_kind = result.gene_kind

            # add to GBK
            gbk_dict[result.gbk_id].genes.append(new_cds)


def translate(feature: SeqFeature, nt_seq: Seq):
    """Translate a CDS sequence if CDS is valid

    Args:
        feature (SeqFeature): CDS SeqFeature
        nt_seq (Seq): CDS nucleotide sequence

    Returns:
        _type_: Seq
    """

    warnings.filterwarnings(action="ignore", category=BiopythonWarning)

    transl_table = "Standard"

    if "transl_table" in feature.qualifiers.keys():
        transl_table = feature.qualifiers["transl_table"][0]

    transl_nt_seq = nt_seq.translate(table=transl_table, to_stop=True)

    return str(transl_nt_seq)


def trim_fuzzy(
    feature: SeqFeature, nt_seq: Seq, fuzzy_start: bool, fuzzy_end: bool, remainder: int
):
    """Method to trim an out-of-frame sequence based on fuzzies

    Fuzzy start: assume that end of sequence is in frame
    and trim from the beginning

    Fuzzy end: assume that the start of sequence is in frame
    and trim from the end

    Args:
        feature (SeqFeature): CDS SeqFeature
        nt_seq (Seq): CDS nucleotide sequence

    Returns:
        _type_: Seq or None
    """

    if fuzzy_start == fuzzy_end:
        raise ValueError(
            "Trimming not possible if fuzzy start and end have the same value"
        )

    if fuzzy_start:
        if remainder == 1:
            trimmed_nt_seq = nt_seq[1:]
        else:
            trimmed_nt_seq = nt_seq[2:]

    if fuzzy_end:
        if remainder == 1:
            trimmed_nt_seq = nt_seq[:-1]
        else:
            trimmed_nt_seq = nt_seq[:-2]

    return trimmed_nt_seq


def check_translation(aa_seq: Seq, nt_seq: Seq, feature: SeqFeature):
    """Check if antiSMASH provided translation matches the one
    generated by Biopython

    Args:
        aa_seq (Seq): antiSMASH translation
        nt_seq (Seq): CDS nucleotide sequence
        feature (SeqFeature): CDS SeqFeature

    Returns:
        _type_: Bool
    """

    transl_nt_seq = translate(feature, nt_seq)

    if str(aa_seq) == str(transl_nt_seq):
        return True

    # case where starting codon may have been set to Methionine
    if str(aa_seq)[1:] == str(transl_nt_seq)[1:]:
        return True

    return False


def get_translation(feature: SeqFeature, nt_seq: Seq):
    """Method to generate a Biopython translation from a CDS sequence

    If nt_seq is not divisible by 3, we try to trim it (trimming is possible
    if either fuzzy start or fuzzy end are present, but not if both or
    neither is present)

    Furthermore, no checks are performed for the presence of start/stop
    codons (biopython cds=True), as there is no influence on downstream analysis


    Args:
        feature (SeqFeature): CDS SeqFeature
        nt_seq (Seq): CDS nucleotide sequence

    Returns:
        _type_: Seq or None
    """

    remainder = len(nt_seq) % 3

    if remainder == 0:
        # if len nt_seq divisible by 3, we translate direclty
        transl_nt_seq = translate(feature, nt_seq)
        return transl_nt_seq

    # if len nt_seq not divisible by 3, we try to trim
    fuzzy_start = False
    fuzzy_end = False
    if str(feature.location.start)[0] in "<>":
        fuzzy_start = True

    if str(feature.location.end)[0] in "<>":
        fuzzy_end = True

    if fuzzy_start and fuzzy_end:
        # if both start and end are fuzzy, trimming will not generate
        # an in-frame sequence
        return None

    if not fuzzy_start and not fuzzy_end:
        # this should not happen, but just in case
        # if no fuzzies, then can't trim using the method below
        return None

    trimmed_nt_seq = trim_fuzzy(feature, nt_seq, fuzzy_start, fuzzy_end, remainder)

    transl_nt_seq = translate(feature, trimmed_nt_seq)

    return transl_nt_seq
