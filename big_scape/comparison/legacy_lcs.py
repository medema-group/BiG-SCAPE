"""Contains code to perform legacy LCS detection as implemented in BiG-SCAPE 1.0"""

# from python
from typing import Any
from difflib import Match, SequenceMatcher

# from other modules
from big_scape.genbank import CDS


def legacy_find_lcs(
    list_a: list[Any], list_b: list[Any]
) -> tuple[Match, list[Match]]:  # pragma no cover
    """Detect longest common substring using sequencematcher

    Args:
        list_a (list[Any]): A list of hashable objects
        list_b (list[Any]): Second list of hashable objects

    Returns:
        tuple[int, int, int]: start, stop and length of longest common substring
    """
    seqmatch = SequenceMatcher(None, list_a, list_b)
    match = seqmatch.find_longest_match(0, len(list_a), 0, len(list_b))
    matching_blocks = seqmatch.get_matching_blocks()
    return match, matching_blocks


def legacy_find_cds_lcs(
    a_cds: list[CDS], b_cds: list[CDS]
) -> tuple[int, int, int, int, bool]:  # pragma no cover
    """Find the longest stretch of matching domains between two CDS lists

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    # forward
    match, matching_blocks = legacy_find_lcs(a_cds, b_cds)
    a_start_fwd = match[0]
    b_start_fwd = match[1]
    fwd_match_len = match[2]

    # reverse. matching blocks is not used in 1.0 logic
    match, matching_blocks_rev = legacy_find_lcs(a_cds, b_cds[::-1])
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
        a_stop = a_start_rev + rev_match_len

        # previously we flipped the entire B array in order to find the LCS
        # now we want to go back to a start and stop position for the unflipped
        # array
        b_start = b_start_rev
        b_stop = b_start + rev_match_len

    # from here on they're the same length
    elif a_cds[a_start_fwd].strand == b_cds[b_start_fwd].strand:
        reverse = False
        a_start = a_start_fwd
        a_stop = a_start_fwd + fwd_match_len

        b_start = b_start_fwd
        b_stop = b_start_fwd + fwd_match_len

    # case where slice length is 1. use the one with the most domains from the seq matc
    elif fwd_match_len == 1:
        max_domains = 0
        for a_idx, b_idx, match_len in matching_blocks:
            if match_len == 0:
                break

            a_domain_count = len(a_cds[a_idx].hsps)
            if a_domain_count <= max_domains:
                continue

            max_domains = a_domain_count
            a_start = a_idx
            a_stop = a_start + 1

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

    return a_start, a_stop, b_start, b_stop, reverse
