"""Contains class to describe domain hits and alignments"""

# from python
from __future__ import annotations
from typing import Optional, TYPE_CHECKING

import tqdm

# from other modules
from big_scape.data import DB


# circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import CDS  # imported in genbank.gbk
    from big_scape.genbank.gbk import GBK


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

        if len(accession) > 10:
            raise ValueError(f"Accession {accession} is not in the expected format")

        return accession

    @staticmethod
    def load_all(gbk_list: list[GBK]) -> None:
        """Load all HSPs and HSP alignments from the database

        This function adds the HSP objects to the CDS objects in cds_list

        Args:
            cds_list (list[CDS]): list of CDS objects that were previously loaded from
            to populate with HSP objects from the database
        """
        # cds_dict = {cds._db_id: cds for cds in cds_list}

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        cds_table = DB.metadata.tables["cds"]
        hsp_table = DB.metadata.tables["hsp"]
        hsp_alignment_table = DB.metadata.tables["hsp_alignment"]

        cds_query = cds_table.select().add_columns(
            cds_table.c.id, cds_table.c.gbk_id, cds_table.c.orf_num
        )

        cursor_result = DB.execute(cds_query)

        cds_orf_to_id = {}
        cds_id_to_hsp: dict[int, list] = {}

        for result in cursor_result.all():
            cds_orf_to_id[(result.gbk_id, result.orf_num)] = result.id
            cds_id_to_hsp[result.id] = []

        hsp_query = (
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
            .join(
                hsp_alignment_table,
                hsp_alignment_table.c.hsp_id == hsp_table.c.id,
            )
        )

        cursor_result = DB.execute(hsp_query)

        for result in cursor_result.all():
            cds_id_to_hsp[result.cds_id].append(result)

        # hsp_alignment_query = hsp_alignment_table.select().add_columns(
        #     hsp_alignment_table.c.hsp_id,
        #     hsp_alignment_table.c.alignment,
        # )

        # cursor_result = DB.execute(hsp_alignment_query)

        # hsp_id_to_alignment_id = {}

        # for result in cursor_result.all():
        #     hsp_id_to_alignment_id[result.hsp_id] = result.id

        progress = tqdm.tqdm(
            gbk_list, desc="Adding db ids to GBK gene data", unit="GBK"
        )

        for gbk in progress:
            for cds in gbk.genes:
                cds._db_id = cds_orf_to_id[(cds.parent_gbk._db_id, cds.orf_num)]

                for hsp_result in cds_id_to_hsp[cds._db_id]:
                    new_hsp = HSP(
                        cds,
                        hsp_result.accession,
                        hsp_result.bit_score,
                        hsp_result.env_start,
                        hsp_result.env_stop,
                    )
                    new_hsp._db_id = hsp_result.id
                    new_hsp.alignment = HSPAlignment(new_hsp, hsp_result.alignment)

                    cds.hsps.append(new_hsp)

        # hsp_select_query = (
        #     hsp_table.select()
        #     .add_columns(
        #         hsp_table.c.id,
        #         hsp_table.c.cds_id,
        #         hsp_table.c.accession,
        #         hsp_table.c.env_start,
        #         hsp_table.c.env_stop,
        #         hsp_table.c.bit_score,
        #         hsp_alignment_table.c.alignment,
        #     )
        #     .join(hsp_alignment_table, hsp_alignment_table.c.hsp_id == hsp_table.c.id)
        #     .where(hsp_table.c.cds_id.in_(cds_dict))
        #     .compile()
        # )

        # cursor_result = DB.execute(hsp_select_query)

        # for result in cursor_result.all():
        #     cds = cds_dict[result.cds_id]
        #     new_hsp = HSP(
        #         cds,
        #         result.accession,
        #         result.bit_score,
        #         result.env_start,
        #         result.env_stop,
        #     )
        #     new_hsp.alignment = HSPAlignment(new_hsp, result.alignment)

        #     cds.hsps.append(new_hsp)

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
