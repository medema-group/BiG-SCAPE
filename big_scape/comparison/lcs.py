"""Contains functions for computing the longest common subsequence between two lists
of domains or CDS in a RecordPair"""

# from python
from difflib import Match, SequenceMatcher

# from other modules
import big_scape.genbank as bs_genbank


def find_lcs(list_a: list, list_b: list) -> tuple[Match, list[Match]]:
    """Detect longest common substring using sequencematcher

    Args:
        list_a (list[T]): A list of hashable objects
        list_b (list[T]): Second list of hashable objects

    Returns:
        tuple[int, int, int]: start, stop and length of longest common substring
    """
    seqmatch = SequenceMatcher(None, list_a, list_b)
    match = seqmatch.find_longest_match(0, len(list_a), 0, len(list_b))
    matching_blocks = seqmatch.get_matching_blocks()
    return match, matching_blocks


def find_cds_lcs(
    a_cds: list[bs_genbank.CDS], b_cds: list[bs_genbank.CDS]
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two CDS lists

    Args:
        pair (RecordPair): A RecordPair object

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    # forward
    match, matching_blocks = find_lcs(a_cds, b_cds)
    a_start_fwd = match[0]
    b_start_fwd = match[1]
    fwd_match_len = match[2]

    # reverse. matching blocks is not used in 1.0 logic
    match, matching_blocks_rev = find_lcs(a_cds, b_cds[::-1])
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

        return a_start, a_stop, b_start, b_stop, reverse

    if rev_larger:
        reverse = True
        a_start = a_start_rev
        a_stop = a_start_rev + rev_match_len

        b_start = len(b_cds) - b_start_rev - rev_match_len
        b_stop = len(b_cds) - b_start_rev

        return a_start, a_stop, b_start, b_stop, reverse

    # equal lengths

    # length of 1. use the matching block with the most domains
    # if all equal, this returns the first
    # TODO: should probably return most central
    if fwd_match_len == 1:
        max = 0
        for a_idx, b_idx, match_len in matching_blocks:
            num_domains = len(a_cds[a_idx : a_idx + match_len])  # noqa
            if num_domains > max:
                max = num_domains

                reverse = False

                a_start = a_idx
                a_stop = a_idx + match_len

                b_start = b_idx
                b_stop = b_idx + match_len

        for a_idx, b_idx, match_len in matching_blocks_rev:
            num_domains = len(a_cds[a_idx : a_idx + match_len])  # noqa
            if num_domains > max:
                max = num_domains

                reverse = True

                a_start = a_idx
                a_stop = a_idx + match_len

                b_start = len(b_cds) - b_idx - match_len
                b_stop = len(b_cds) - b_idx

        return a_start, a_stop, b_start, b_stop, reverse

    # equal length, but not 1
    # default to forward
    # default to first match
    # TODO: should probably return most central
    reverse = False
    a_start = a_start_fwd
    a_stop = a_start_fwd + fwd_match_len

    b_start = b_start_fwd
    b_stop = b_start_fwd + fwd_match_len

    return a_start, a_stop, b_start, b_stop, reverse
