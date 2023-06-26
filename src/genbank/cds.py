""" Module containing code to load and store antiSMASH CDSs"""

# from python
from __future__ import annotations
import logging
from typing import Optional, TYPE_CHECKING

# from dependencies
from Bio.SeqFeature import SeqFeature
from sortedcontainers import SortedList

# from other modules
from src.errors import InvalidGBKError
from src.data import DB
from src.hmm import HSP


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import GBK  # imported earlier in src.file_input.load_files


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
        hsps: SortedList[HSP]
    """

    def __init__(self, nt_start: int, nt_stop: int):
        self.nt_start = nt_start
        self.nt_stop = nt_stop
        self.orf_num: Optional[int] = None
        self.parent_gbk: Optional[GBK] = None
        self.gene_kind: Optional[str] = None
        self.strand: Optional[int] = None
        self.aa_seq: str = ""
        self.hsps: SortedList[HSP] = SortedList()

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
            hsp (src.hmmer.hsp): The HSP to be added to this CDS
            overlap_cutoff (float, optional): cutoff threshold for overlap. Defaults to
            0.1

        """
        # if no hsps added yet, just add and continue
        if len(self.hsps) == 0:
            self.hsps.add(new_hsp)
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
        self.hsps.add(new_hsp)

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
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None):
        """Creates a cds object from a region feature in a GBK file

        Note that this will not add biosynthetic information if no parent gbk is passed
        """

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

        # add parent if it exists
        if parent_gbk is None:
            return cds

        cds.parent_gbk = parent_gbk
        if parent_gbk.as_version == "4":
            if "sec_met" not in feature.qualifiers:
                cds.gene_kind = ""
                return cds

            for sec_met_value in feature.qualifiers["sec_met"]:
                if "Kind" in sec_met_value:
                    cds.gene_kind = sec_met_value[6:]  # trim "Kind: "
                    return cds

        if "gene_kind" in feature.qualifiers:
            cds.gene_kind = str(feature.qualifiers["gene_kind"][0])

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
