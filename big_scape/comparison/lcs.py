"""Contains functions for computing the longest common subsequence between two lists
of domains or CDS in a RecordPair"""

# from python
from typing import Any
from difflib import Match, SequenceMatcher

# from other modules
import big_scape.genbank as bs_genbank


def find_lcs(list_a: list[Any], list_b: list[Any]) -> tuple[Match, list[Match]]:
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


def find_protocore_distance(protocluster: bs_genbank.Region, idx: int) -> int:
    """Find the distance between a CDS and the closest protocore

    Args:
        protocluster (Region): protocluster
        idx (int): index of CDS

    Returns:
        int: distance to closest protocore
    """

    if not isinstance(protocluster, bs_genbank.ProtoCluster):
        raise TypeError("protocluster must be a protocluster")

    min_dist = None
    for protocore_idx in protocluster.proto_core_cds_idx:
        dist = abs(protocore_idx - idx)
        if min_dist is None or dist < min_dist:
            min_dist = dist
    return min_dist


def find_cds_lcs_region(
    a_cds: list[bs_genbank.CDS], b_cds: list[bs_genbank.CDS]
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two CDS lists

    If there are LCS of the same length, the LCS closest to the middle of the region
    is preferred (TODO)
    TODO: maybe this is not useful at all

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    # forward
    match, matching_blocks = find_lcs(a_cds, b_cds)
    a_start_fwd = match[0]
    b_start_fwd = match[1]
    fwd_match_len = match[2]

    # reverse
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


def find_domain_lcs_region(
    a_cds: list[bs_genbank.CDS], b_cds: list[bs_genbank.CDS]
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two lists of domains

    This takes CDS as arguments, but uses the domains within the CDS to find the LCS

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """

    a_domains = []
    b_domains = []
    for cds in a_cds:
        a_domains.extend(cds.hsps)
    for cds in b_cds:
        b_domains.extend(cds.hsps)
    # forward
    match, matching_blocks = find_lcs(a_domains, b_domains)
    a_start_fwd = match[0]
    b_start_fwd = match[1]
    fwd_match_len = match[2]

    # reverse
    match, matching_blocks_rev = find_lcs(a_domains, b_domains[::-1])
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

        b_start = len(b_domains) - b_start_rev - rev_match_len
        b_stop = len(b_domains) - b_start_rev

        return a_start, a_stop, b_start, b_stop, reverse

    # equal lengths
    # match of length 1 means we pick something in the middle
    if fwd_match_len == 1:
        # first find which region is shorter in terms of cds
        a_cds_len = len(a_cds)
        b_cds_len = len(b_cds)

        # default to A. if B is shorter, use B
        if a_cds_len <= b_cds_len:
            use_cds = a_cds
            use_domains = a_domains
            matching_block_idx = 0
        else:
            use_cds = b_cds
            use_domains = b_domains
            matching_block_idx = 1

        # generate a CDS to index dict
        cds_idx_dict = {cds: i for i, cds in enumerate(use_cds)}

        # go through all LCS matches and find the one with the most central CDS
        middle = len(use_cds) / 2
        min = None
        for matching_block in matching_blocks:
            # I don't even know why there is a match of len 0 when there are matches
            # of len 1
            if matching_block[2] == 0:
                continue

            idx = matching_block[matching_block_idx]

            domain = use_domains[idx]
            cds_idx = cds_idx_dict[domain.cds]

            # find the distance to the middle
            distance = abs(middle - cds_idx)

            if min is None or distance < min:
                min = distance
                a_start = matching_block[0]
                a_stop = matching_block[0] + matching_block[2]
                b_start = matching_block[1]
                b_stop = matching_block[1] + matching_block[2]

        return a_start, a_stop, b_start, b_stop, False

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
