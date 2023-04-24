"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from __future__ import annotations
from typing import Type
from difflib import SequenceMatcher

# from this module
from .binning import BGCPair


class LCS:
    """Class describing the longest common subsequence of domains between a pair of BGCs

    Properties:
        a_start: int
        b_start: int
        match_len: int
        reverse: bool
    """

    def __init__(self, a_start: int, b_start: int, match_len: int, reverse: bool):
        self.a_start = a_start
        self.b_start = b_start
        self.match_len = match_len
        self.reverse = reverse

    @classmethod
    def get_pair_domain_lcs(cls: Type[LCS], pair: BGCPair) -> LCS:
        """Retrieve the longest common subsequence of domains for a pair of BGC records
        Also returns whether the sequence is reversed

        Args:
            pair (BGCPair): Pair of BGCs to retrieve the LCS for

        Returns:
            tuple[list[str], bool]: List of strings that represents the LCS of domains and
            a boolean indicating whether the LCS is was found in the reverse sequence
        """
        a_hsps = pair.region_a.get_hsps()
        b_hsps = pair.region_b.get_hsps()

        a_domains = [a_hsp.domain for a_hsp in a_hsps]
        b_domains = [b_hsp.domain for b_hsp in b_hsps]

        # forward
        seqmatch = SequenceMatcher(None, a_domains, b_domains)
        match = seqmatch.find_longest_match(0, len(a_domains), 0, len(b_domains))
        a_start_fwd = match[0]
        b_start_fwd = match[1]
        fwd_match_len = match[2]

        # reverse
        seqmatch = SequenceMatcher(None, a_domains, b_domains[::-1])
        match = seqmatch.find_longest_match(0, len(a_domains), 0, len(b_domains))
        a_start_rev = match[0]
        b_start_rev = match[1]
        rev_match_len = match[2]

        if fwd_match_len >= rev_match_len:
            reverse = False
            a_start = a_start_fwd
            b_start = b_start_fwd
            match_len = fwd_match_len
        else:
            reverse = True
            a_start = a_start_rev
            b_start = b_start_rev
            match_len = rev_match_len

        return cls(a_start, b_start, match_len, reverse)

    def __eq__(self, __value: object) -> bool:
        """Checks whether this LCS object and another object are equal

        Args:
            __value (object): Object to compare this LCS object to

        Returns:
            bool: True if this LCS object is equal to another (LCS) object
        """
        if not isinstance(__value, LCS):
            return False

        conditions = [
            self.a_start == __value.a_start,
            self.b_start == __value.b_start,
            self.match_len == __value.match_len,
            self.reverse == __value.reverse,
        ]
        return all(conditions)

    def __repr__(self) -> str:
        """Returns a human-readable string representation of an instance of this class

        Returns:
            str: Human-readable class description
        """
        return (
            f"LCS: as: {self.a_start}, bs: {self.b_start}, len: {self.match_len} "
            f"reverse: {self.reverse}"
        )
