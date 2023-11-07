"""Contains functions for computing the longest common subsequence between two lists
of domains or CDS in a RecordPair

Into any of the main functions in this module goes a RecordPair object, which contains
two regions. These regions can be either protoclusters or full regions. The functions
in this module are used to find the longest common subsequence between the two regions
in the RecordPair. This is used in extension as a "seed" for the extension.

The result from the main functions are tuples with the following structure:
    (a_start, a_stop, b_start, b_stop, reverse)
Where a and b are the two regions in the RecordPair. Start is inclusive, stop is
exclusive. Reverse is a boolean indicating whether the match is in reverse
"""

# from python
import logging
from typing import Any
from difflib import Match, SequenceMatcher

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.comparison as bs_comparison
import big_scape.hmm as bs_hmm


def find_lcs(list_a: list[Any], list_b: list[Any]) -> tuple[Match, list[Match]]:
    """Detect longest common substring using sequencematcher

    Args:
        list_a (list[T]): A list of hashable objects
        list_b (list[T]): Second list of hashable objects

    Returns:
        tuple[int, int, int]: start, stop and length of longest common substring
    """
    seqmatch = SequenceMatcher(None, list_a, list_b, False)
    match = seqmatch.find_longest_match(0, len(list_a), 0, len(list_b))
    matching_blocks = seqmatch.get_matching_blocks()
    return match, matching_blocks


def find_protocore_distance(region: bs_genbank.ProtoCluster, idx: int) -> int:
    """Find the distance between a CDS and the closest protocore

    Args:
        protocluster (Region): protocluster
        idx (int): index of CDS

    Returns:
        int: distance to closest protocore
    """

    if not isinstance(region, bs_genbank.ProtoCluster):
        raise TypeError("region must be a protocluster")

    min_dist = None
    for protocore_idx in region.proto_core_cds_idx:
        dist = abs(protocore_idx - idx)
        if min_dist is None or dist < min_dist:
            min_dist = dist

    if min_dist is None:
        raise ValueError("No protocore found")  # pragma: no cover # should never happen

    return min_dist


def get_lcs_protocores(
    pair: bs_comparison.RecordPair, matching_blocks: list[Match], reverse: bool
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two protocluster regions,
    preferring matches which are closest to a protocore

    Returns the cds indexes of start and stop, whether the match is in reverse,
    and whether the match is in the protocore

    Args:
        pair (RecordPair): RecordPair object
        matching_blocks (list[Match]): list of matching blocks
        reverse (bool): whether the match is in reverse

    Returns:
        tuple[int, int, int, int, bool, bool]: a_start, a_stop, b_start, b_stop,
        match in protocore
    """

    if not isinstance(pair.region_a, bs_genbank.ProtoCluster):
        raise TypeError("region_a must be a protocluster")

    if not isinstance(pair.region_b, bs_genbank.ProtoCluster):
        raise TypeError("region_b must be a protocluster")

    a_min_dist = None
    b_min_dist = None
    a_best = None
    b_best = None

    # we need to know if the domains are in the protocore, but all we have is the
    # pair.region_b.proto_core_cds_idx. we need to make a similar index on domain level
    a_proto_core_domain_idx = set()
    b_proto_core_domain_idx = set()

    domain_idx = 0
    for idx, cds in enumerate(pair.region_a.get_cds()):
        if idx not in pair.region_a.proto_core_cds_idx:
            domain_idx += len(cds.hsps)
            continue

        for _ in cds.hsps:
            a_proto_core_domain_idx.add(domain_idx)
            domain_idx += 1

    domain_idx = 0
    for idx, cds in enumerate(pair.region_b.get_cds()):
        if idx not in pair.region_b.proto_core_cds_idx:
            domain_idx += len(cds.hsps)
            continue

        for _ in cds.hsps:
            b_proto_core_domain_idx.add(domain_idx)
            domain_idx += 1

    # we need to keep track of this to correct for the reverse later
    b_num_domains = domain_idx

    for a_idx, b_idx, match_len in matching_blocks:
        # if match len > 1, check all the indexes in the match

        if reverse:
            b_idx = b_num_domains - b_idx - match_len

        if match_len > 1:
            a_in_protocore = any(
                [
                    idx in a_proto_core_domain_idx
                    for idx in range(a_idx, a_idx + match_len)
                ]
            )
            b_in_protocore = any(
                [
                    idx in b_proto_core_domain_idx
                    for idx in range(b_idx, b_idx + match_len)
                ]
            )
        else:
            a_in_protocore = a_idx in a_proto_core_domain_idx
            b_in_protocore = b_idx in b_proto_core_domain_idx

        # exit early if both are in a protocore
        if a_in_protocore and b_in_protocore:
            # flip b_idx again
            if reverse:
                b_idx = b_num_domains - b_idx - match_len
            return a_idx, a_idx + match_len, b_idx, b_idx + match_len, True

        # from this point we can assume we need to find the distance to the closest
        # protocore

        # if match_len > 1, use whichever is closest to a protocore
        if match_len > 1:
            left_dist = find_protocore_distance(pair.region_a, a_idx)
            right_dist = find_protocore_distance(pair.region_a, a_idx + match_len - 1)
            a_dist = min(left_dist, right_dist)

            left_dist = find_protocore_distance(pair.region_b, b_idx)
            right_dist = find_protocore_distance(pair.region_b, b_idx + match_len - 1)
            b_dist = min(left_dist, right_dist)
        else:
            a_dist = find_protocore_distance(pair.region_a, a_idx)
            b_dist = find_protocore_distance(pair.region_b, b_idx)

        a_better = False
        if a_min_dist is None or a_dist < a_min_dist:
            a_better = True

        b_better = False
        if b_min_dist is None or b_dist < b_min_dist:
            b_better = True

        if a_better and b_better:
            a_min_dist = a_dist
            b_min_dist = b_dist
            a_best = a_idx
            b_best = b_idx
            best_len = match_len

    # now we have the best pair of indexes, or the first one which is fine

    if a_best is None or b_best is None:
        raise ValueError("No match found")

    a_start = a_best
    a_stop = a_best + best_len
    b_start = b_best
    b_stop = b_best + best_len

    return a_start, a_stop, b_start, b_stop, False


def find_domain_lcs_region(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two lists of domains

    This takes CDS as arguments, but uses the domains within the CDS to find the LCS

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    logging.debug("region lcs")

    # these are regions, so we can get the full range of CDS
    a_cds = pair.region_a.get_cds(True)
    b_cds = pair.region_b.get_cds(True)

    # get lists of domains and assemble a dictionary of domain idx to cds idx
    # so that we can return the cds indexes later
    a_domains: list[bs_hmm.HSP] = []
    a_domain_cds_idx = {}

    b_domains: list[bs_hmm.HSP] = []
    b_domain_cds_idx = {}

    for cds_idx, cds in enumerate(a_cds):
        for i in range(len(a_domains), len(a_domains) + len(cds.hsps)):
            a_domain_cds_idx[i] = cds_idx
        a_domains.extend(cds.hsps)

    for cds_idx, cds in enumerate(b_cds):
        for i in range(len(b_domains), len(b_domains) + len(cds.hsps)):
            b_domain_cds_idx[i] = cds_idx
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

    if fwd_match_len == 0 and rev_match_len == 0:
        logging.error(
            "No match found in LCS. This should not happen after first jaccard"
        )
        logging.error("a domains: %s", a_domains)
        logging.error("b domains: %s", b_domains)
        raise RuntimeError("No match found in LCS.")

    if fwd_larger:
        reverse = False
        a_start = a_start_fwd
        a_stop = a_start_fwd + fwd_match_len

        b_start = b_start_fwd
        b_stop = b_start_fwd + fwd_match_len

        a_cds_start = a_domain_cds_idx[a_start]
        a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
        b_cds_start = b_domain_cds_idx[b_start]
        b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

        return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, False

    if rev_larger:
        reverse = True
        a_start = a_start_rev
        a_stop = a_start_rev + rev_match_len

        b_start = len(b_domains) - b_start_rev - rev_match_len
        b_stop = len(b_domains) - b_start_rev

        a_cds_start = a_domain_cds_idx[a_start]
        a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
        b_cds_start = b_domain_cds_idx[b_start]
        b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

        return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, False

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

        a_cds_start = a_domain_cds_idx[a_start]
        a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
        b_cds_start = b_domain_cds_idx[b_start]
        b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

        return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, False

    # equal length, but not 1
    # default to forward
    # default to first match
    # TODO: should probably return most central

    reverse = False
    a_start = a_start_fwd
    a_stop = a_start_fwd + fwd_match_len

    b_start = b_start_fwd
    b_stop = b_start_fwd + fwd_match_len

    a_cds_start = a_domain_cds_idx[a_start]
    a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
    b_cds_start = b_domain_cds_idx[b_start]
    b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

    return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse


def find_domain_lcs_protocluster(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two protocluster records,
    using domains

    Args:
        pair (RecordPair): RecordPair object

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    logging.debug("pc lcs")

    # we really need protoclusters here
    if not isinstance(pair.region_a, bs_genbank.ProtoCluster):
        raise TypeError("region_a must be a protocluster")

    if not isinstance(pair.region_b, bs_genbank.ProtoCluster):
        raise TypeError("region_b must be a protocluster")

    a_cds = pair.region_a.get_cds()
    b_cds = pair.region_b.get_cds()

    a_domains = []
    b_domains = []
    for idx, cds in enumerate(a_cds):
        a_domains.extend(cds.hsps)
    for idx, cds in enumerate(b_cds):
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

    # forward
    forward_lcs = get_lcs_protocores(pair, matching_blocks, False)
    in_protocore_fwd = forward_lcs[4]

    # reverse
    reverse_lcs = get_lcs_protocores(pair, matching_blocks_rev, True)
    in_protocore_rev = reverse_lcs[4]

    # if a match is found both in reverse and forward that contains protocores in both
    # regions, use the longest match. if matches are equal length, use forward
    if in_protocore_fwd and in_protocore_rev:
        a_start_fwd, a_stop_fwd = forward_lcs[0:2]
        a_start_rev, a_stop_rev = reverse_lcs[0:2]
        if a_stop_fwd - a_start_fwd >= a_stop_rev - a_start_rev:
            reverse = False
            return forward_lcs[0:4] + (reverse,)
        else:
            reverse = True
            return reverse_lcs[0:4] + (reverse,)

    # if a match is found in forward, use that
    if in_protocore_fwd:
        reverse = False
        return forward_lcs[0:4] + (reverse,)

    # if a match is found in reverse, use that
    if in_protocore_rev:
        reverse = True
        return reverse_lcs[0:4] + (reverse,)

    # if no match is found in either, use the longest match
    if fwd_match_len >= rev_match_len:
        reverse = False
        return (
            a_start_fwd,
            a_start_fwd + fwd_match_len,
            b_start_fwd,
            b_start_fwd + fwd_match_len,
            reverse,
        )
    else:
        reverse = True
        return (
            a_start_rev,
            a_start_rev + rev_match_len,
            b_start_rev,
            b_start_rev + rev_match_len,
            reverse,
        )
