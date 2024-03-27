"""Contains class to describe domain hits and alignments"""

# from python
from __future__ import annotations
from typing import Optional, TYPE_CHECKING

# from other modules
from big_scape.data import DB


# circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import CDS  # imported in genbank.gbk


class HSP:
    """Describes a CDS - Domain relationship"""

    def __init__(
        self, cds: CDS, domain: str, score: float, env_start: int, env_stop: int
    ) -> None:
        self.cds = cds
        self.domain = HSP.sanitize_accession(domain)
        self.score = score
        self.env_start = env_start
        self.env_stop = env_stop

        self.alignment: Optional[HSPAlignment] = None

        # db specific fields
        self._db_id: Optional[int] = None

    def save(self, commit=True) -> None:
        """Saves this hsp object to a database and optionally executes a commit

        Args:
            commit (bool, optional): Whether to commit immediately after inserting this
            CDS. Defaults to True."""

        parent_cds_id = None
        if self.cds is not None and self.cds._db_id is not None:
            parent_cds_id = self.cds._db_id

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        hsp_table = DB.metadata.tables["hsp"]
        insert_query = (
            hsp_table.insert()
            .returning(hsp_table.c.id)
            .values(
                cds_id=parent_cds_id,
                accession=self.domain,
                env_start=self.env_start,
                env_stop=self.env_stop,
                bit_score=self.score,
            )
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

    def __repr__(self) -> str:
        return f"{self.domain} ({self.score:.2f})"

    def __gt__(self, __o: object) -> bool:
        if not isinstance(__o, HSP):
            raise NotImplementedError()

        if self.cds.nt_start < __o.cds.nt_start:
            return False

        if self.env_start < __o.env_start:
            return False
        if self.env_start > __o.env_start:
            return True

        # TODO: sorting is possibly based on e-value, not on score
        # not that it should matter that much, but domain order affects AI
        if self.score < __o.score:
            return True

        return False

    def __eq__(self, __o: object) -> bool:
        """Return whether this HSP object and another object are equal

        Will always return false if the object compared is not an instance of HSP or str

        If compared to a string, checks if this HSP accession is equal to that string

        Args:
            __o (object): Target object to compare to

        Returns:
            bool: True if __o and self are equal
        """
        # special case if we are comparing this to a string
        if isinstance(__o, str):
            # we do not care about version numbers in this comparison, so strip it
            return self.domain == __o

        if not isinstance(__o, HSP):
            raise NotImplementedError()

        return __o.domain == self.domain

    def __hash__(self) -> int:
        return hash(self.domain)

    @staticmethod
    def sanitize_accession(accession: str):
        if len(accession) < 8:
            return accession

        if accession[7] == ".":
            return accession[:7]

        if len(accession) > 9:
            raise ValueError(f"Accession {accession} is not in the expected format")

        return accession

    @staticmethod
    def load_all(cds_list: list[CDS]) -> None:
        """Load all HSPs and HSP alignments from the database

        This function adds the HSP objects to the CDS objects in cds_list

        Args:
            cds_list (list[CDS]): list of CDS objects that were previously loaded from
            to populate with HSP objects from the database
        """
        cds_dict = {cds._db_id: cds for cds in cds_list}

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        hsp_table = DB.metadata.tables["hsp"]
        hsp_alignment_table = DB.metadata.tables["hsp_alignment"]

        hsp_select_query = (
            hsp_table.select()
            .add_columns(
                hsp_table.c.id,
                hsp_table.c.cds_id,
                hsp_table.c.accession,
                hsp_table.c.env_start,
                hsp_table.c.env_stop,
                hsp_table.c.bit_score,
                hsp_alignment_table.c.alignment,
            )
            .join(hsp_alignment_table, hsp_alignment_table.c.hsp_id == hsp_table.c.id)
            .where(hsp_table.c.cds_id.in_(cds_dict))
            .compile()
        )

        cursor_result = DB.execute(hsp_select_query)

        for result in cursor_result.all():
            cds = cds_dict[result.cds_id]
            new_hsp = HSP(
                cds,
                result.accession,
                result.bit_score,
                result.env_start,
                result.env_stop,
            )
            new_hsp.alignment = HSPAlignment(new_hsp, result.alignment)

            cds.hsps.append(new_hsp)

    @staticmethod
    def has_overlap(hsp_a: HSP, hsp_b: HSP) -> bool:
        """Return True if there is overlap between this

        Args:
            cds_b (CDS): CDS to compare

        Returns:
            bool: whether there is overlap between this cds and another
        """
        has_overlap = HSP.len_overlap(hsp_a, hsp_b) > 0
        return has_overlap

    @staticmethod
    def len_overlap(hsp_a: HSP, hsp_b: HSP) -> int:
        """Return the length of the overlap between this CDS and another

        Args:
            cds_b (CDS): CDS to compare

        Returns:
            int: length of the overlap between this CDS and another
        """

        if hsp_a.env_start < hsp_b.env_start:
            left = hsp_b.env_start
        else:
            left = hsp_a.env_start

        if hsp_a.env_stop > hsp_b.env_stop:
            right = hsp_b.env_stop
        else:
            right = hsp_a.env_stop

        # limit to > 0
        return max(0, right - left)

    @staticmethod
    def overlap_filter(hsp_list: list[HSP], overlap_cutoff=0.1) -> list[HSP]:
        """Filters overlapping HSP

        Args:
            hsp_list (list[HSP]): a list of high scoring protein hits
            overlap_cutoff: The maximum percentage sequence length domains can overlap
                without being discarded

        Returns:
            list[HSP]: a filtered list of high scoring protein hits
        """
        # for this the hsps have to be sorted by any of the start coordinates
        # this is for the rare case where two hsps have the same bit score.
        # to replicate BiG-SCAPE output, the hsp with a lower start coordinate
        # will end up being chosen
        sorted_hsps = sorted(hsp_list, key=lambda hsp: hsp.cds.nt_start)

        filtered_list = []
        for i in range(len(sorted_hsps) - 1):
            for j in range(i + 1, len(sorted_hsps)):
                hsp_a: HSP = sorted_hsps[i]
                hsp_b: HSP = sorted_hsps[j]

                # check if we are the same CDS
                if hsp_a == hsp_b:
                    # check if there is overlap between the domains
                    if not HSP.has_overlap(hsp_a, hsp_b):
                        continue

                    len_aa_overlap = HSP.len_overlap(hsp_a, hsp_b)

                    # calculate percentage of total hsp length is the aa overlap between
                    # two hsp
                    overlap_perc_a = len_aa_overlap / (
                        hsp_a.cds.nt_stop - hsp_a.cds.nt_start
                    )
                    overlap_perc_b = len_aa_overlap / (
                        hsp_b.cds.nt_stop - hsp_b.cds.nt_start
                    )
                    # check if the amount of overlap is significant
                    if (
                        overlap_perc_a < overlap_cutoff
                        and overlap_perc_b < overlap_cutoff
                    ):
                        continue
                    # rounding with 1 decimal to mirror domtable files
                    hsp_a_score = round(float(hsp_a.score), 1)
                    hsp_b_score = round(float(hsp_b.score), 1)
                    if hsp_a_score >= hsp_b_score:  # see which has a better score
                        filtered_list.append(hsp_a)
                    elif hsp_a_score < hsp_b_score:
                        filtered_list.append(hsp_b)

        return filtered_list


class HSPAlignment:
    """Describes a HSP - Domain alignment"""

    def __init__(self, hsp: HSP, align_string: str) -> None:
        self.hsp = hsp
        self.align_string = align_string

        # database specific fields
        self._db_id: Optional[int] = None

    def save(self, commit=True) -> None:
        """Saves this hsp alignment object to a database and optionally executes a
        commit

        Args:
            commit (bool, optional): Whether to commit immediately after inserting this
            CDS. Defaults to True."""

        parent_hsp_id = None
        if self.hsp is not None and self.hsp._db_id is not None:
            parent_hsp_id = self.hsp._db_id

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        hsp_align_table = DB.metadata.tables["hsp_alignment"]
        insert_query = hsp_align_table.insert().values(
            hsp_id=parent_hsp_id,
            alignment=self.align_string,
        )

        DB.execute(insert_query, commit)

    def __repr__(self) -> str:
        return (
            f"Alignment of hsp {str(self.hsp)} and domain {self.hsp.domain}: "
            "{self.alignment}"
        )

    def __eq__(self, __o: object) -> bool:
        """Return whether this HSP alignment and another object are equal

        Will always return false if the object compared is not an instance of
        HSPAlignment

        Args:
            __o (object): Target object to compare to

        Returns:
            bool: True if __o and self are equal
        """
        if not isinstance(__o, HSPAlignment):
            return False

        conditions = [
            self.hsp == __o.hsp,
            self.align_string == __o.align_string,
        ]
        return all(conditions)
