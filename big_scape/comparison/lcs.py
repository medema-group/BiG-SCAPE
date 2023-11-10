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

NOTE: The matches correspond to slices of region CDS that do not have domains!
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

    NOTE: The LCS correspond to slices of region CDS that do not have domains!
    These slices are converted to full CDS ranges after extend, not here in LCS

    Args:
        pair (RecordPair): RecordPair object
        matching_blocks (list[Match]): list of matching blocks
        reverse (bool): whether the match is in reverse

    Returns:
        tuple[int, int, int, int, bool, bool]: a_start, a_stop, b_start, b_stop,
        match in protocore
    """

    if not isinstance(pair.record_a, bs_genbank.ProtoCluster):
        raise TypeError("record_a must be a protocluster")

    if not isinstance(pair.record_b, bs_genbank.ProtoCluster):
        raise TypeError("record_b must be a protocluster")

    a_min_dist = None
    b_min_dist = None
    a_best = None
    b_best = None

    # we need to know if the domains are in the protocore, but all we have is the
    # pair.record_b.proto_core_cds_idx. we need to make a similar index on domain level
    a_proto_core_domain_idx = set()
    b_proto_core_domain_idx = set()

    domain_idx = 0
    for idx, cds in enumerate(pair.record_a.get_cds_with_domains()):
        if idx not in pair.record_a.proto_core_cds_idx:
            domain_idx += len(cds.hsps)
            continue

        for _ in cds.hsps:
            a_proto_core_domain_idx.add(domain_idx)
            domain_idx += 1

    domain_idx = 0
    for idx, cds in enumerate(pair.record_b.get_cds_with_domains()):
        if idx not in pair.record_b.proto_core_cds_idx:
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
            left_dist = find_protocore_distance(pair.record_a, a_idx)
            right_dist = find_protocore_distance(pair.record_a, a_idx + match_len - 1)
            a_dist = min(left_dist, right_dist)

            left_dist = find_protocore_distance(pair.record_b, b_idx)
            right_dist = find_protocore_distance(pair.record_b, b_idx + match_len - 1)
            b_dist = min(left_dist, right_dist)
        else:
            a_dist = find_protocore_distance(pair.record_a, a_idx)
            b_dist = find_protocore_distance(pair.record_b, b_idx)

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


def find_bio_or_middle_lcs(
    a_cds: list[bs_genbank.CDS],
    b_cds: list[bs_genbank.CDS],
    a_domains: list[bs_hmm.HSP],
    b_domains: list[bs_hmm.HSP],
    matching_blocks: list[Match],
    matching_blocks_rev: list[Match],
    a_domain_cds_idx: dict[int, int],
    b_domain_cds_idx: dict[int, int],
) -> tuple[int, int, int, int, bool]:
    """Find the most central match out of all LCS matches, or a match containing a
    biosynthetic gene

    This is done by first selecting the shorter record in terms of CDS, and then
    finding the match that is closest to the middle of the CDS

    Args:
        a_cds (list[CDS]): List of CDS for A
        b_cds (list[CDS]): List of CDS for B
        a_domains (list[HSP]): List of domains for A
        b_domains (list[HSP]): List of domains for B
        matching_blocks (list[Match]): List of matching blocks
        matching_blocks_rev (list[Match]): List of matching blocks in reverse
        a_domain_cds_idx (dict[int, int]): Dictionary of domain idx to cds idx for A
        b_domain_cds_idx (dict[int, int]): Dictionary of domain idx to cds idx for B

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """

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
            reverse = False

    # go through all reverse, too
    for matching_block in matching_blocks_rev:
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
            b_start = len(b_domains) - matching_block[1] - matching_block[2]
            b_stop = len(b_domains) - matching_block[1]
            reverse = True

    a_cds_start = a_domain_cds_idx[a_start]
    a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
    b_cds_start = b_domain_cds_idx[b_start]
    b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

    return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse


def find_domain_lcs_region(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two lists of domains

    This takes CDS as arguments, but uses the domains within the CDS to find the LCS

    NOTE: The LCS correspond to slices of region CDS that do not have domains!
    These slices need to be converted to full CDS ranges later

    Approach:
    - Pick the largest LCS with a biosynthetic gene
    - If there are multiple LCS with a biosynthetic gene of equal length, pick the one
        that is closest to the middle of the region
    - If there are no LCS with a biosynthetic gene, pick the largest
    - If there are multiple LCS of equal length, pick the one that is closest to the
        middle of the region

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    logging.debug("region lcs")

    # these are regions, so we can get the full range of CDS
    a_cds = pair.record_a.get_cds_with_domains(True)
    b_cds = pair.record_b.get_cds_with_domains(True)

    # working on domains, not cds
    a_domains: list[bs_hmm.HSP] = []
    b_domains: list[bs_hmm.HSP] = []

    # dictionary of domain index to cds index to quickly find the cds of a domain
    a_domain_cds_idx = {}
    b_domain_cds_idx = {}

    # list of domain idx whose genes are biosynthetic, to quickly find if a domain is
    # part of a biosynthetic gene
    a_biosynthetic_domain_cds: list[int] = []
    b_biosynthetic_domain_cds: list[int] = []

    for cds_idx, cds in enumerate(a_cds):
        for i in range(len(a_domains), len(a_domains) + len(cds.hsps)):
            a_domain_cds_idx[i] = cds_idx
        a_domains.extend(cds.hsps)
        if cds.gene_kind == "biosynthetic":
            a_biosynthetic_domain_cds.append(cds_idx)

    for cds_idx, cds in enumerate(b_cds):
        for i in range(len(b_domains), len(b_domains) + len(cds.hsps)):
            b_domain_cds_idx[i] = cds_idx
        b_domains.extend(cds.hsps)
        if cds.gene_kind == "biosynthetic":
            b_biosynthetic_domain_cds.append(cds_idx)

    # forward
    match, matching_blocks_fwd = find_lcs(a_domains, b_domains)
    fwd_match_len = match[2]

    # reverse
    match, matching_blocks_rev = find_lcs(a_domains, b_domains[::-1])
    rev_match_len = match[2]

    # quickly check if we didn't find an LCS
    if fwd_match_len == 0 and rev_match_len == 0:
        logging.error(
            "No match found in LCS. This should not happen after first jaccard"
        )
        logging.error("a domains: %s", a_domains)
        logging.error("b domains: %s", b_domains)
        raise RuntimeError("No match found in LCS.")

    # now we need to do something silly. we want to assemble a list of these matching
    # blocks, but we want to keep track of whether they are in reverse or not. this will
    # make it so we can do everything we need to do in one loop later
    matchin_block_dirs = []
    for matching_block in matching_blocks_fwd:
        matchin_block_dirs.append((matching_block + (False,)))

    for matching_block in matching_blocks_rev:
        matchin_block_dirs.append((matching_block + (True,)))

    # this is where the fun begins. we will use these lists to decide which match to
    # return later

    # tuple is idx, length, reverse
    longest_biosynthetic: list[tuple[int, int, bool]] = []
    longest: list[tuple[int, int, bool]] = []
    # tuple is idx, distance to middle, reverse
    central_biosynthetic: list[tuple[int, int, bool]] = []
    central: list[tuple[int, int, bool]] = []

    for match_idx, matching_block_dir in enumerate(matchin_block_dirs):
        start_a = matching_block_dir[0]
        stop_a = matching_block_dir[0] + matching_block_dir[2]
        start_b = matching_block_dir[1]
        stop_b = matching_block_dir[1] + matching_block_dir[2]
        length = matching_block_dir[2]
        reverse = matching_block_dir[3]

        # I don't understand why, but zero-length blocks exist sometimes. skip them
        if length == 0:
            continue

        # fix b start and stop if in reverse
        if reverse:
            start_b = len(b_domains) - start_b - length
            stop_b = len(b_domains) - stop_b

        # check if the match contains a biosynthetic gene
        has_biosynthetic = False
        for biosynthetic_idx in a_biosynthetic_domain_cds:
            if start_a <= biosynthetic_idx < stop_a:
                has_biosynthetic = True
                break

        for biosynthetic_idx in b_biosynthetic_domain_cds:
            if start_b <= biosynthetic_idx < stop_b:
                has_biosynthetic = True
                break

        # select for biosynthetic or normal list
        use_longest_list = longest_biosynthetic if has_biosynthetic else longest
        use_central_list = central_biosynthetic if has_biosynthetic else central

        # Length

        # clear the list if it's not empty and the current match is longer
        if len(use_longest_list) > 0 and length > use_longest_list[0][1]:
            use_longest_list.clear()

        # add the match to the list if it's empty or the current match is equal length
        # or longer than the existing match
        if len(use_longest_list) == 0 or length >= use_longest_list[0][1]:
            use_longest_list.append((match_idx, length, reverse))

        # distance to middle

        # use the shorter cds list to determine the distance to middle
        use_cds = a_cds if len(a_cds) <= len(b_cds) else b_cds
        use_idx = a_domain_cds_idx if len(a_cds) <= len(b_cds) else b_domain_cds_idx
        use_start = start_a if len(a_cds) <= len(b_cds) else start_b
        use_stop = stop_a if len(a_cds) <= len(b_cds) else stop_b

        middle = len(use_cds) / 2
        # calculate the distance from either side of the match to the middle
        distance = round(
            min(
                abs(middle - use_idx[use_start]),
                abs(middle - use_idx[use_stop - 1]),
            )
        )

        # clear the list if it's not empty and the current match is closer to the middle
        if len(use_central_list) > 0 and distance < use_central_list[0][1]:
            use_central_list.clear()

        # add the match to the list if it's empty or the current match is equal distance
        # or closer to the middle than the existing match
        if len(use_central_list) == 0 or distance <= use_central_list[0][1]:
            use_central_list.append((match_idx, distance, reverse))

    # now we have everything we need. we need to decide which match to return
    # just go top to bottom and return the first match in the list
    # remember that the lists are [match_idx, length/dist, reverse]
    # refer to docstring for decision making here
    if len(longest_biosynthetic) == 1:
        match_idx = longest_biosynthetic[0][0]

    elif len(central_biosynthetic) > 0:
        match_idx = central_biosynthetic[0][0]

    elif len(longest) == 1:
        match_idx = longest[0][0]

    elif len(central) > 0:
        match_idx = central[0][0]

    else:
        # this should never happen
        raise RuntimeError("No match found in LCS.")

    relevant_match = matchin_block_dirs[match_idx]
    a_start = relevant_match[0]
    a_stop = relevant_match[0] + relevant_match[2]
    b_start = relevant_match[1]
    b_stop = relevant_match[1] + relevant_match[2]
    reverse = relevant_match[3]

    # fix b start and stop if in reverse
    if reverse:
        old_start = b_start
        b_start = len(b_domains) - b_stop
        b_stop = len(b_domains) - old_start

    a_cds_start = a_domain_cds_idx[a_start]
    a_cds_stop = a_domain_cds_idx[a_stop]
    b_cds_start = b_domain_cds_idx[b_start]
    b_cds_stop = b_domain_cds_idx[b_stop]

    return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse


def find_domain_lcs_protocluster(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two protocluster records,
    using domains

    NOTE: The LCS correspond to slices of region CDS that do not have domains!
    These slices need to be converted to full CDS ranges later

    Args:
        pair (RecordPair): RecordPair object

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """
    logging.debug("pc lcs")

    # we really need protoclusters here
    if not isinstance(pair.record_a, bs_genbank.ProtoCluster):
        raise TypeError("record_a must be a protocluster")

    if not isinstance(pair.record_b, bs_genbank.ProtoCluster):
        raise TypeError("record_b must be a protocluster")

    a_cds = pair.record_a.get_cds_with_domains()
    b_cds = pair.record_b.get_cds_with_domains()

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
