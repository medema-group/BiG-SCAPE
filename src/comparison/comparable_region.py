"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from __future__ import annotations
import logging
from typing import Type, TYPE_CHECKING, Optional
from difflib import SequenceMatcher

# from dependencies

# from other modules
from src.genbank import BGCRecord


# from this module

# from other modules
from src.hmm import HSP


# from circular imports
if TYPE_CHECKING:
    from .binning import BGCPair  # not sure why this one throws a circular import error


class ComparableRegion:
    """Class describing the comparable region between a pair of BGCs

    Properties:
        a_start: int
        b_start: int
        a_len: int
        b_len: int
        reverse: bool
        domain_strings: list[str]
    """

    def __init__(
        self,
        pair: BGCPair,
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

    def get_domain_sets(self, regenerate=False) -> tuple[set[HSP], set[HSP]]:
        """Returns a tuple containing sets of domains within the comparable region of
        two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one

        Returns:
            tuple[set[HSP], set[HSP]]
        """
        if regenerate or self.domain_sets is None:
            a_domain_list, b_domain_list = self.get_domain_lists()

            self.domain_sets = (set(a_domain_list), set(b_domain_list))

        return self.domain_sets

    def get_domain_lists(self, regenerate=False) -> tuple[list[HSP], list[HSP]]:
        """Returns a tuple corresponding (ordered) lists of CDS domains within the
        comparable region of two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one

        Returns:
            tuple[list[HSP], list[HSP]]
        """

        if regenerate or self.domain_lists is None:
            a_start = self.a_start
            a_stop = self.a_stop

            b_start = self.b_start
            b_stop = self.b_stop

            a_cds_list = self.pair.region_a.get_cds()[a_start:a_stop]
            b_cds_list = self.pair.region_b.get_cds()[b_start:b_stop]

            if self.reverse:
                b_cds_list = b_cds_list[::-1]

            a_domain_list = []
            for a_cds in a_cds_list:
                a_domain_list.extend(a_cds.hsps)

            b_domain_list = []
            for b_cds in b_cds_list:
                b_domain_list.extend(b_cds.hsps)

            self.domain_lists = (a_domain_list, b_domain_list)

        return self.domain_lists

    @classmethod
    def create_domain_lcs(
        cls: Type[ComparableRegion], pair: BGCPair
    ) -> ComparableRegion:
        """Retrieve the longest common subsequence of domains for a pair of BGC records
        Also returns whether the sequence is reversed

        Args:
            pair (BGCPair): Pair of BGCs to retrieve the LCS for

        Returns:
            tuple[list[str], bool]: List of strings that represents the LCS of domains and
            a boolean indicating whether the LCS is was found in the reverse sequence
        """
        # the idea here is that each cds has a list of domains that are matched against
        # we concatenate the domains within a CDS, and the list of concatenated domains
        # for all cds is called the domain string
        a_cds = pair.region_a.get_cds()
        b_cds = pair.region_b.get_cds()

        # forward
        seqmatch = SequenceMatcher(None, a_cds, b_cds)
        match = seqmatch.find_longest_match(0, len(a_cds), 0, len(b_cds))
        a_start_fwd = match[0]
        b_start_fwd = match[1]
        fwd_match_len = match[2]

        # reverse
        seqmatch = SequenceMatcher(None, a_cds, b_cds[::-1])
        match = seqmatch.find_longest_match(0, len(a_cds), 0, len(b_cds))
        a_start_rev = match[0]
        b_start_rev = match[1]
        rev_match_len = match[2]

        if fwd_match_len >= rev_match_len:
            reverse = False
            a_start = a_start_fwd
            a_stop = a_start_fwd + fwd_match_len

            b_start = b_start_fwd
            b_stop = b_start_fwd + fwd_match_len
        else:
            reverse = True
            a_start = a_start_rev
            a_stop = a_start_rev + rev_match_len

            # previously we flipped the entire B array in order to find the LCS
            # now we want to go back to a start and stop position for the unflipped
            # array
            b_start = b_start_rev
            b_stop = b_start + rev_match_len

        return cls(pair, a_start, a_stop, b_start, b_stop, reverse)

    def log_comparable_region(self, label="<") -> None:
        """Prints a debug level log of the comparable region

        Args:
            label (str): A string to visually indicate the comparable region
        """
        if logging.getLogger().level > logging.DEBUG:
            return

        a_cds_list = self.pair.region_a.get_cds()
        b_cds_list = self.pair.region_b.get_cds()

        if self.reverse:
            b_cds_list = b_cds_list[::-1]

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

            log_line = f"{i:<3} {a_domains:<25} {a_region_str:<10} {b_domains:<25} {b_region_str:<10}"
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
        record: BGCRecord, cds_start: int, cds_stop: int, end_inclusive=False
    ) -> bool:
        stop = cds_stop
        if end_inclusive:
            stop += 1

        for cds in record.get_cds(True)[cds_start:stop]:
            if "biosynthetic" in cds.gene_kind:
                return True

        return False
