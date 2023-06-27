"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from __future__ import annotations
import logging
from typing import TYPE_CHECKING, Optional
from difflib import SequenceMatcher

# from other modules
from src.genbank import BGCRecord


# from this module

# from other modules
from src.hmm import HSP


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
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
        self.domain_dicts: Optional[
            tuple[dict[HSP, list[int]], dict[HSP, list[int]]]
        ] = None

    def get_domain_sets(
        self, regenerate=False, cache=True
    ) -> tuple[set[HSP], set[HSP]]:
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
        regenerate is set to True or Cache is set to false.

        If reverse is set to none, uses the reverse propoerty on comparable region

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one,
            cache (bool): whether to cache the result for faster retrieval
            reverse (bool, optional): Whether to return a reversed list for region B. Defaults to
            None.

        Returns:
            tuple[SortedList[HSP], SortedList[HSP]]

        """

        if regenerate or self.domain_lists is None:
            a_start = self.a_start
            a_stop = self.a_stop

            b_start = self.b_start
            b_stop = self.b_stop

            a_cds_list = self.pair.region_a.get_cds()[a_start:a_stop]
            b_cds_list = self.pair.region_b.get_cds(reverse=self.reverse)[
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

    def find_lcs(self) -> tuple[int, int, int, int]:
        """Retrieve the longest common subsequence of domains for a pair of BGC records

        Args:
            pair (BGCPair): Pair of BGCs to retrieve the LCS for

        Returns:
            tuple[int]: tuple of integers corresponding to:
            [a_lcs_start, a_lcs_stop, b_lcs_start, b_lcs_stop]
        """
        # the idea here is that each cds has a list of domains that are matched against
        # we concatenate the domains within a CDS, and the list of concatenated domains
        # for all cds is called the domain string
        a_cds = self.pair.region_a.get_cds()
        b_cds = self.pair.region_b.get_cds()

        # forward
        seqmatch = SequenceMatcher(None, a_cds, b_cds)
        match = seqmatch.find_longest_match(0, len(a_cds), 0, len(b_cds))
        a_start_fwd = match[0]
        b_start_fwd = match[1]
        fwd_match_len = match[2]

        # reverse
        rev_seqmatch = SequenceMatcher(None, a_cds, b_cds[::-1])
        match = rev_seqmatch.find_longest_match(0, len(a_cds), 0, len(b_cds))
        a_start_rev = match[0]
        b_start_rev = match[1]
        rev_match_len = match[2]

        fwd_larger = fwd_match_len > rev_match_len
        rev_larger = fwd_match_len < rev_match_len

        if fwd_larger:
            reverse = False
            a_start = a_start_fwd
            a_stop = a_start_fwd + fwd_match_len

            b_start = b_start_fwd
            b_stop = b_start_fwd + fwd_match_len

        elif rev_larger:
            reverse = True
            a_start = a_start_rev
            a_stop = a_start_rev + rev_match_len + 1

            # previously we flipped the entire B array in order to find the LCS
            # now we want to go back to a start and stop position for the unflipped
            # array
            b_start = b_start_rev
            b_stop = b_start + rev_match_len + 1

        # from here on they're the same length
        elif (
            self.pair.region_a.get_cds()[a_start_fwd].strand
            == self.pair.region_b.get_cds()[b_start_fwd].strand
        ):
            reverse = False
            a_start = a_start_fwd
            a_stop = a_start_fwd + fwd_match_len

            b_start = b_start_fwd
            b_stop = b_start_fwd + fwd_match_len

        # case where slice length is 1. use the one with the most domains from the seq matc
        elif fwd_match_len == 1:
            max_domains = 0
            for a_idx, b_idx, match_len in seqmatch.get_matching_blocks():
                if match_len == 0:
                    break

                a_domain_count = len(a_cds[a_idx].hsps)
                if a_domain_count <= max_domains:
                    continue

                max_domains = a_domain_count
                a_start = a_idx
                a_stop = a_start

                # TODO: only the else should trigger at this point
                if a_cds[a_idx].strand == b_cds[b_idx].strand:
                    b_start = b_idx
                    b_stop = b_start
                    reverse = False
                else:
                    b_start = len(b_cds) - b_idx - 1
                    b_stop = b_start + 1
                    reverse = True

        # default to taking forward lcs
        else:
            a_start = a_start_fwd
            a_stop = a_start_fwd + fwd_match_len
            b_start = b_start_fwd
            b_stop = b_start + fwd_match_len
            reverse = False

        self.a_start = a_start
        self.a_stop = a_stop
        self.b_start = b_start
        self.b_stop = b_stop
        self.reverse = reverse

        return (a_start, a_stop, b_start, b_stop)

    def log_comparable_region(self, label="<") -> None:  # pragma: no cover
        """Prints a debug level log of the comparable region

        Args:
            label (str): A string to visually indicate the comparable region
        """
        if logging.getLogger().level > logging.DEBUG:
            return

        a_cds_list = self.pair.region_a.get_cds()
        b_cds_list = self.pair.region_b.get_cds(reverse=self.reverse)

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
        record: BGCRecord,
        cds_start: int,
        cds_stop: int,
        end_inclusive=False,
        reverse=False,
    ) -> bool:
        stop = cds_stop
        if end_inclusive:
            stop += 1

        for cds in record.get_cds(True, reverse=reverse)[cds_start:stop]:
            if cds.gene_kind is None:
                continue

            if "biosynthetic" in cds.gene_kind:
                return True

        return False
