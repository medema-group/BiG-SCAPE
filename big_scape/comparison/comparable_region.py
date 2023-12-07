"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from __future__ import annotations
from typing import TYPE_CHECKING

# from other modules
from big_scape.genbank import BGCRecord

# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from .record_pair import (
        RecordPair,
    )  # not sure why this one throws a circular import error


class ComparableRegion:
    """Class describing the comparable region between a pair of BGCs

    Important notes:

    This class is used to store the comparable region between two records. It is used to
    calculate the JC, AI, and DSS scores between two records. The start and stop of the
    comparable region correspond to the CDS indexes of each record.

    Comparable regions are initialized with coordinates CDS including domains. This is used in LCS
    and extend to determine the comparable region. After LCS and extend, the comparable region is
    recalculated to include all CDS, including those without domains using the inflate method.

    Reverse indicates whether, before working on the CDS lists, the CDS list for record B should be
    reversed. This is detected in LCS. The B CDS start and stop do not change, and reflect the start
    and stop once the list is already reversed.

    Properties:
        pair: BGCPair
        a_start: int
        a_start: int
        b_start: int
        b_stop: int
        domain_a_start: int
        domain_b_start: int
        domain_a_stop: int
        domain_b_stop: int
        reverse: bool
    """

    def __init__(
        self,
        a_start: int,
        a_stop: int,
        b_start: int,
        b_stop: int,
        domain_a_start: int,
        domain_b_start: int,
        domain_a_stop: int,
        domain_b_stop: int,
        reverse: bool,
    ):
        # store lcs cds without any extensions
        self.lcs_a_start = a_start
        self.lcs_b_start = b_start
        self.lcs_a_stop = a_stop
        self.lcs_b_stop = b_stop

        # store lcs domain
        self.lcs_domain_a_start = domain_a_start
        self.lcs_domain_b_start = domain_b_start
        self.lcs_domain_a_stop = domain_a_stop
        self.lcs_domain_b_stop = domain_b_stop

        # store possibly extended comparable region cds
        self.a_start = a_start
        self.b_start = b_start
        self.a_stop = a_stop
        self.b_stop = b_stop

        # store possibly extended comparable region domain
        self.domain_a_start = domain_a_start
        self.domain_b_start = domain_b_start
        self.domain_a_stop = domain_a_stop
        self.domain_b_stop = domain_b_stop

        self.reverse = reverse

    def inflate_a(self, record_a: BGCRecord) -> None:
        """Inflates the A region start and stop to include all CDS, including those
        without domains.

        This should be done once after LCS/Extend. There is no way to know if this
        comparable region has already been inflated, so this method should only be
        called once.
        """
        a_cds_list = record_a.get_cds()
        a_cds_with_domains = record_a.get_cds_with_domains()
        lcs_a_start_orf_num = a_cds_with_domains[self.lcs_a_start].orf_num
        lcs_a_stop_orf_num = a_cds_with_domains[self.lcs_a_stop - 1].orf_num
        a_start_orf_num = a_cds_with_domains[self.a_start].orf_num
        a_stop_orf_num = a_cds_with_domains[self.a_stop - 1].orf_num

        for idx, cds in enumerate(a_cds_list):
            if cds.orf_num == lcs_a_start_orf_num:
                self.lcs_a_start = idx

            if cds.orf_num == a_start_orf_num:
                self.a_start = idx

            if cds.orf_num == lcs_a_stop_orf_num:
                self.lcs_a_stop = idx + 1

            if cds.orf_num == a_stop_orf_num:
                self.a_stop = idx + 1
                break

    def inflate_b(self, record_b: BGCRecord) -> None:
        """Inflates the B region start and stop to include all CDS, including those
        without domains.

        This should be done once after LCS/Extend. There is no way to know if this
        comparable region has already been inflated, so this method should only be
        called once.
        """
        b_cds_list = record_b.get_cds()
        b_cds_with_domains = record_b.get_cds_with_domains()

        if self.reverse:
            b_cds_list = b_cds_list[::-1]
            b_cds_with_domains = b_cds_with_domains[::-1]

        lcs_b_start_orf_num = b_cds_with_domains[self.lcs_b_start].orf_num
        lcs_b_stop_orf_num = b_cds_with_domains[self.lcs_b_stop - 1].orf_num
        b_start_orf_num = b_cds_with_domains[self.b_start].orf_num
        b_stop_orf_num = b_cds_with_domains[self.b_stop - 1].orf_num

        for idx, cds in enumerate(b_cds_list):
            if cds.orf_num == lcs_b_start_orf_num:
                self.lcs_b_start = idx

            if cds.orf_num == b_start_orf_num:
                self.b_start = idx

            if cds.orf_num == lcs_b_stop_orf_num:
                self.lcs_b_stop = idx + 1

            if cds.orf_num == b_stop_orf_num:
                self.b_stop = idx + 1
                break

    def inflate(self, pair: RecordPair):
        """Inflates the comparable region to include all CDS, including those without
        domains.

        This should be done once after LCS/Extend. There is no way to know if this
        comparable region has already been inflated, so this method should only be
        called once.
        """
        self.inflate_a(pair.record_a)
        self.inflate_b(pair.record_b)

    def to_tuple(self):
        """Returns a tuple representation of this comparable region

        Returns:
            tuple: Tuple representation of this comparable region
        """
        return (
            self.lcs_a_start,
            self.lcs_a_stop,
            self.lcs_b_start,
            self.lcs_b_stop,
            self.a_start,
            self.a_stop,
            self.b_start,
            self.b_stop,
            self.reverse,
            self.lcs_domain_a_start,
            self.lcs_domain_a_stop,
            self.lcs_domain_b_start,
            self.lcs_domain_b_stop,
        )

    def __eq__(self, __o: object) -> bool:
        """Checks whether this LCS object and another object are equal

        Args:
            __value (object): Object to compare this LCS object to

        Returns:
            bool: True if this LCS object is equal to another (LCS) object
        """
        if not isinstance(__o, ComparableRegion):
            raise NotImplementedError()

        conditions = [
            self.a_start == __o.a_start,
            self.b_start == __o.b_start,
            self.a_stop == __o.a_stop,
            self.a_stop == __o.a_stop,
            self.reverse == __o.reverse,
        ]
        return all(conditions)

    def __repr__(self) -> str:
        """Returns a human-readable string representation of an instance of this class

        Returns:
            str: Human-readable class description
        """
        if self.reverse:
            reverse_string = "B is reversed"
        else:
            reverse_string = "B is not reversed"
        b_start = self.b_start
        b_stop = self.b_stop

        return (
            f"Comparable region: A {self.a_start}-{self.a_stop}, "
            f"B {b_start}-{b_stop}, {reverse_string}"
        )

    @staticmethod
    def cds_range_contains_biosynthetic(
        record: BGCRecord,
        cds_start: int,
        cds_stop: int,
        end_inclusive=False,
        reverse=False,
    ) -> bool:
        """Return true if a range of cds within this record contains a gene that is
        marked as biosynthetic

        Args:
            record (BGCRecord): BGC record object
            cds_start (int): which cds index to start searching in
            cds_stop (int): Which cds index to stop searching
            end_inclusive (bool, optional): Whether to include the cds at cds_stop in
            the search. Defaults to False.
            reverse (bool, optional): Whether to reverse the list before applying
            cds_start and cds_stop. Defaults to False.

        Returns:
            bool: True if the given range contains a biosynthetic gene
        """
        stop = cds_stop
        if end_inclusive:
            stop += 1

        for cds in record.get_cds_with_domains(reverse=reverse)[cds_start:stop]:
            if cds.gene_kind is None:
                continue

            if "biosynthetic" in cds.gene_kind:
                return True

        return False
