"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Optional

# from other modules
from big_scape.genbank import BGCRecord
from big_scape.hmm import HSP

# from this module
from .legacy_lcs import legacy_find_cds_lcs


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from .binning import (
        RecordPair,
    )  # not sure why this one throws a circular import error


class ComparableRegion:
    """Class describing the comparable region between a pair of BGCs

    Properties:
        pair: BGCPair
        a_start: int
        a_start: int
        b_start: int
        b_stop: int
        reverse: bool
        domain_lists: Optional[tuple[list[HSP], list[HSP]]]
        domain_sets: Optional[tuple[set[HSP], set[HSP]]]
        domain_dicts: Optional[tuple[dict[HSP, list[int]], dict[HSP, list[int]]]]
    """

    def __init__(
        self,
        pair: RecordPair,
        a_start: int,
        a_stop: int,
        b_start: int,
        b_stop: int,
        reverse: bool,
    ):
        self.pair = pair
        self.a_start = a_start
        self.b_start = b_start

        self.a_stop = a_stop
        self.b_stop = b_stop

        self.reverse = reverse

        self.domain_lists: Optional[tuple[list[HSP], list[HSP]]] = None
        self.domain_sets: Optional[tuple[set[HSP], set[HSP]]] = None
        self.domain_dicts: Optional[
            tuple[dict[HSP, list[int]], dict[HSP, list[int]]]
        ] = None

    def get_domain_sets(
        self, regenerate=False, cache=True
    ) -> tuple[set[HSP], set[HSP]]:
        """Returns a tuple containing sets of domains within the comparable region of
        two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True. This is done to minimize the computational cost of
        intializing a set from a list at the cost of memory.

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one
            cache (bool): whether to cache the result of this operation

        Returns:
            tuple[set[HSP], set[HSP]]
        """

        if regenerate or self.domain_sets is None:
            a_domain_list, b_domain_list = self.get_domain_lists(cache=cache)

            if not cache:
                return (set(a_domain_list), set(b_domain_list))

            self.domain_sets = (set(a_domain_list), set(b_domain_list))

        return self.domain_sets

    def get_domain_lists(
        self, regenerate=False, cache=True, reverse=None
    ) -> tuple[list[HSP], list[HSP]]:
        """Returns a tuple corresponding (ordered) lists of CDS domains within the
        comparable region of two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True. This is done to minimize the computational cost of
        intializing a set from a list at the cost of memory.

        If reverse is set to none, uses the reverse propoerty on comparable region

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one,
            cache (bool): whether to cache the result for faster retrieval
            reverse (bool, optional): Whether to return a reversed list for region B.
            Defaults to None.

        Returns:
            tuple[list[HSP], list[HSP]]
        """

        if regenerate or self.domain_lists is None:
            a_start = self.a_start
            a_stop = self.a_stop

            b_start = self.b_start
            b_stop = self.b_stop

            a_cds_list = self.pair.region_a.get_cds_with_domains()[a_start:a_stop]
            b_cds_list = self.pair.region_b.get_cds_with_domains(reverse=self.reverse)[
                b_start:b_stop
            ]

            a_domain_list: list[HSP] = []
            for a_cds in a_cds_list:
                a_domain_list.extend(a_cds.hsps)

            b_domain_list: list[HSP] = []
            for b_cds in b_cds_list:
                b_domain_list.extend(b_cds.hsps)

            if not cache:
                return (a_domain_list, b_domain_list)

            self.domain_lists = (a_domain_list, b_domain_list)

        return self.domain_lists

    def get_domain_dicts(
        self, regenerate=False
    ) -> tuple[dict[HSP, list[int]], dict[HSP, list[int]]]:
        """Returns a dictionary of domains for each BGC in this comaprable region.
        Dictionary keys are domains and ints are the index of that domain in the domain
        list from which the dictionary was generated. This will always be the list
        returned by self.get_domain_list at the moment this method is called

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one

        Returns:
            tuple[dict[HSP, list[int]], dict[HSP, list[int]]]: dictionary of domain
            accessions to list indexes
        """

        if regenerate or self.domain_dicts is None:
            domain_dict_a: dict[HSP, list[int]] = {}
            domain_dict_b: dict[HSP, list[int]] = {}

            domain_list_a, domain_list_b = self.get_domain_lists()

            for idx_a, domain_a in enumerate(domain_list_a):
                if domain_a not in domain_dict_a:
                    domain_dict_a[domain_a] = []

                domain_dict_a[domain_a].append(idx_a)

            for idx_b, domain_b in enumerate(domain_list_b):
                if domain_b not in domain_dict_b:
                    domain_dict_b[domain_b] = []

                domain_dict_b[domain_b].append(idx_b)

            self.domain_dicts = (domain_dict_a, domain_dict_b)

        return self.domain_dicts

    def find_lcs(self):
        """Retrieve the longest common subsequence of domains for a pair of BGC records

        Args:
            pair (BGCPair): Pair of BGCs to retrieve the LCS for

        Returns:
            tuple[int]: tuple of integers corresponding to:
            [a_lcs_start, a_lcs_stop, b_lcs_start, b_lcs_stop]
        """
        a_cds = self.pair.region_a.get_cds_with_domains()
        b_cds = self.pair.region_b.get_cds_with_domains()

        a_start, a_stop, b_start, b_stop, reverse = legacy_find_cds_lcs(a_cds, b_cds)

        self.a_start = a_start
        self.a_stop = a_stop
        self.b_start = b_start
        self.b_stop = b_stop
        self.reverse = reverse

    def log_comparable_region(self, label="<") -> None:  # pragma: no cover
        """Prints a debug level log of the comparable region

        Args:
            label (str): A string to visually indicate the comparable region
        """
        if logging.getLogger().level > logging.DEBUG:
            return

        a_cds_list = self.pair.region_a.get_cds_with_domains()
        b_cds_list = self.pair.region_b.get_cds_with_domains(reverse=self.reverse)

        b_start = self.b_start
        b_stop = self.b_stop

        for i in range(max(len(a_cds_list), len(b_cds_list))):
            a_domains = ""
            if i < len(a_cds_list):
                a_domains = "|".join(
                    [hsp.domain[2:].split(".")[0] for hsp in a_cds_list[i].hsps]
                )

            b_domains = ""
            if i < len(b_cds_list):
                b_domains = "|".join(
                    [hsp.domain[2:].split(".")[0] for hsp in b_cds_list[i].hsps]
                )

            a_region_str = ""
            b_region_str = ""

            a_in_region = i > self.a_start and i < self.a_stop
            if a_in_region:
                a_region_str = label
            if i == self.a_start:
                a_region_str = "START"
            if i == self.a_stop:
                a_region_str = "STOP"

            b_in_region = i > b_start and i < b_stop
            if b_in_region:
                b_region_str = label
            if i == b_start:
                b_region_str = "START"
            if i == b_stop:
                b_region_str = "STOP"

            log_line = " ".join(
                [
                    f"{i:<3}",
                    f"{a_domains:<25}",
                    f"{a_region_str:<10}",
                    f"{b_domains:<25}",
                    f"{b_region_str:<10}",
                ]
            )
            logging.debug(log_line)

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

        for cds in record.get_cds_with_domains(True, reverse=reverse)[cds_start:stop]:
            if cds.gene_kind is None:
                continue

            if "biosynthetic" in cds.gene_kind:
                return True

        return False
