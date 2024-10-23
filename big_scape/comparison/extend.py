"""Contains methods for extension of comparable regions with record pairs"""

# from python
import math
import logging

# from other modules
import big_scape.hmm as bs_hmm

# from this module
from .comparable_region import ComparableRegion
from .record_pair import RecordPair


def reset(pair: RecordPair) -> None:
    """Resets the extension of a pair's comparable region

    Args:
        pair: The record pair to reset
    """

    pair.comparable_region.a_start = 0
    pair.comparable_region.b_start = 0
    pair.comparable_region.a_stop = len(list(pair.record_a.get_cds_with_domains()))
    pair.comparable_region.b_stop = len(list(pair.record_b.get_cds_with_domains()))
    pair.comparable_region.domain_a_start = 0
    pair.comparable_region.domain_b_start = 0
    pair.comparable_region.domain_a_stop = len(list(pair.record_a.get_hsps()))
    pair.comparable_region.domain_b_stop = len(list(pair.record_b.get_hsps()))

    pair.comparable_region.reverse = False


def len_check(pair: RecordPair, min_len_perc: float) -> bool:
    """Checks if a pair's comparable region is of sufficient length

    Length is checked based on domains, the comparable region should be longer than
    or equal to a min_len_percentage of the shortest record in the pair

    Args:
        pair: The record pair to check
        min_len_perc: The minimum length of the comparable region
    """
    a_len = pair.comparable_region.domain_a_stop - pair.comparable_region.domain_a_start
    a_total_len = len(list(pair.record_a.get_hsps()))

    b_len = pair.comparable_region.domain_b_stop - pair.comparable_region.domain_b_start
    b_total_len = len(list(pair.record_b.get_hsps()))

    if a_len / a_total_len >= min_len_perc or b_len / b_total_len >= min_len_perc:
        return True

    return False


def biosynthetic_check(pair: RecordPair) -> bool:
    """Checks if any element of a pair contains a biosynthetic domain

    returns true if
    - the comparable region contains a biosynthetic gene

    Args:
        pair: The record pair to check
    """

    if ComparableRegion.cds_range_contains_biosynthetic(
        pair.record_a,
        pair.comparable_region.a_start,
        pair.comparable_region.a_stop,
        False,
        False,
    ):
        return True

    if ComparableRegion.cds_range_contains_biosynthetic(
        pair.record_b,
        pair.comparable_region.b_start,
        pair.comparable_region.b_stop,
        False,
        pair.comparable_region.reverse,
    ):
        return True

    return False


def extend_glocal(pair: RecordPair) -> None:
    """Includes extension of shortest upstream/downstream arms in comparable region

    Args:
        pair (RecordPair): record pair to extend
    """
    # upstream
    if pair.comparable_region.domain_a_start <= pair.comparable_region.domain_b_start:
        pair.comparable_region.a_start = 0
        pair.comparable_region.domain_a_start = 0
    else:
        pair.comparable_region.b_start = 0
        pair.comparable_region.domain_b_start = 0

    # downstream
    a_len = len(list(pair.record_a.get_cds_with_domains()))
    b_len = len(list(pair.record_b.get_cds_with_domains()))
    a_domain_len = len(list(pair.record_a.get_hsps()))
    b_domain_len = len(list(pair.record_b.get_hsps()))
    a_arm = a_domain_len - pair.comparable_region.domain_a_stop
    b_arm = b_domain_len - pair.comparable_region.domain_b_stop
    if a_arm <= b_arm:
        pair.comparable_region.a_stop = a_len
        pair.comparable_region.domain_a_stop = a_domain_len
    else:
        pair.comparable_region.b_stop = b_len
        pair.comparable_region.domain_b_stop = b_domain_len


def extend(
    pair: RecordPair,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist_perc: float,
) -> None:
    """Extends a comparable region

    This will extend the included set of cds in a pair based on a scoring
    mechanism. If the pair in the comparable region consists of protoclusters, the
    this will not be limited to the bounds of those protoclusters

    NOTE: The start-stop slices returned from this function correspond to CDS slices,
    without domains only.

    Args:
        pair: The record pair to extend
        match: The score for a match
        mismatch: The score for a mismatch
        gap: The score for a gap
        max_match_dist_perc: The maximum distance of a matching domain before it is
    """

    logging.debug("before extend:")
    logging.debug(pair.comparable_region)

    # get the cds lists
    # TODO: base extend on all domains in case of protoclusters, allow extend beyond
    # protocluster border
    a_domains = list(pair.record_a.get_hsps())
    b_domains = list(pair.record_b.get_hsps())

    a_max_dist = math.floor(len(a_domains) * max_match_dist_perc)
    b_max_dist = math.floor(len(b_domains) * max_match_dist_perc)

    # reverse b if necessary. This might be true after LCS
    if pair.comparable_region.reverse:
        b_domains = b_domains[::-1]

    # we will try the following approach:
    # the shorter cds is the query
    # the longer cds is the target
    # generate a dictionary of domain positions in the target
    # we try to extend the target by finding matching domains in the query
    # when we do not find a match, this counts as a mismatch and we add a penalty
    # when we find a match, add gap penalties for the number of gaps we have to insert
    # if we find a match before the current position, subtract a gap penalty
    # we will only match when the position from the current extension is below a cutoff
    # length

    if len(a_domains) > len(b_domains):
        query_domains = b_domains
        query_start = pair.comparable_region.b_start
        query_domain_start = pair.comparable_region.domain_b_start
        query_stop = pair.comparable_region.b_stop
        query_domain_stop = pair.comparable_region.domain_b_stop
        target_domains = a_domains
        target_start = pair.comparable_region.a_start
        target_domain_start = pair.comparable_region.domain_a_start
        target_stop = pair.comparable_region.a_stop
        target_domain_stop = pair.comparable_region.domain_a_stop
        max_match_dist = a_max_dist
    else:
        query_domains = a_domains
        query_start = pair.comparable_region.a_start
        query_domain_start = pair.comparable_region.domain_a_start
        query_stop = pair.comparable_region.a_stop
        query_domain_stop = pair.comparable_region.domain_a_stop
        target_domains = b_domains
        target_start = pair.comparable_region.b_start
        target_domain_start = pair.comparable_region.domain_b_start
        target_stop = pair.comparable_region.b_stop
        target_domain_stop = pair.comparable_region.domain_b_stop
        max_match_dist = b_max_dist

    # generate an index of domain positions in the target
    # the lists in this index will be sorted by position, in ascending order
    target_index = get_target_indexes(target_domains)
    query_index = get_query_indexes(query_domains)

    if target_domain_stop != len(target_domains) and query_domain_stop != len(
        query_domains
    ):
        query_exp, query_dom_exp, target_exp, target_dom_exp, score = score_extend(
            query_domains,
            query_index,
            query_stop,
            query_domain_stop,
            target_index,
            target_stop,
            target_domain_stop,
            match,
            mismatch,
            gap,
            max_match_dist,
        )

        # set the new start and stop positions
        if len(a_domains) > len(b_domains):
            pair.comparable_region.b_stop += query_exp
            pair.comparable_region.domain_b_stop += query_dom_exp
            pair.comparable_region.a_stop += target_exp
            pair.comparable_region.domain_a_stop += target_dom_exp
        else:
            pair.comparable_region.a_stop += query_exp
            pair.comparable_region.domain_a_stop += query_dom_exp
            pair.comparable_region.b_stop += target_exp
            pair.comparable_region.domain_b_stop += target_dom_exp

    if target_domain_start != 0 and query_domain_start != 0:
        query_exp, query_dom_exp, target_exp, target_dom_exp, score = score_extend_rev(
            query_domains,
            query_index,
            query_start,
            query_domain_start,
            target_index,
            target_start,
            target_domain_start,
            match,
            mismatch,
            gap,
            max_match_dist,
        )

        # extend left
        if len(a_domains) > len(b_domains):
            pair.comparable_region.b_start -= query_exp
            pair.comparable_region.domain_b_start -= query_dom_exp
            pair.comparable_region.a_start -= target_exp
            pair.comparable_region.domain_a_start -= target_dom_exp
        else:
            pair.comparable_region.a_start -= query_exp
            pair.comparable_region.domain_a_start -= query_dom_exp
            pair.comparable_region.b_start -= target_exp
            pair.comparable_region.domain_b_start -= target_dom_exp

    logging.debug("after extend:")
    logging.debug(pair.comparable_region)


def get_target_indexes(
    target_domains: list[bs_hmm.HSP],
) -> dict[str, list[tuple[int, int]]]:
    """Generate an index of domain cds positions in the target

    The dictionary that is returned will have the following structure:
    {
        domain: [(cds_idx, domain_idx), ...]
    }

    All indexes are sorted by position, in ascending order

    Args:
        target (list[bs_hmm.HSP]): The target domain list

    Returns:
        dict[str, list[tuple[int, int]]]: The target index dictionary
    """
    target_index: dict[str, list[tuple[int, int]]] = {}

    curr_cds_num = target_domains[0].cds.orf_num
    cds_idx = 0
    for domain_idx, hsp in enumerate(target_domains):
        if hsp.cds.orf_num != curr_cds_num:
            cds_idx += 1
            curr_cds_num = hsp.cds.orf_num

        if hsp.domain not in target_index:
            target_index[hsp.domain] = []

        target_index[hsp.domain].append((cds_idx, domain_idx))

    return target_index


def get_query_indexes(query_domains: list[bs_hmm.HSP]) -> dict[int, int]:
    """Generate an index linking each domain to their cds index

    Dictionary that is returned will have the following structure:
    {
        domain_idx: cds_idx
    }

    domain_idx and cds_idx represent the order of domains/cds in the whole record

    Args:
        query (list[bs_genbank.CDS]): the query cds list

    Returns:
        dict[int, int]: The query index dictionary
    """
    query_index: dict[int, int] = {}

    curr_cds_num = query_domains[0].cds.orf_num
    cds_idx = 0
    for domain_idx, hsp in enumerate(query_domains):
        if hsp.cds.orf_num != curr_cds_num:
            cds_idx += 1
            curr_cds_num = hsp.cds.orf_num

        query_index[domain_idx] = cds_idx
    return query_index


def score_extend(
    query_domains: list[bs_hmm.HSP],
    query_index: dict[int, int],
    query_start: int,
    query_domain_start: int,
    target_index: dict[str, list[tuple[int, int]]],
    target_start: int,
    target_domain_start: int,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist: int,
) -> tuple[int, int, int, int, int]:
    """
    Calculate the score for extending a query sequence to a target sequence.

    TODO: This is a copy of score_extend_rev. if possible, refactor to remove
    duplication

    Args:
        query (list): A list of Cds objects representing the full query record cds
        query_index (dict): A dictionary linking each domain index to its cds index
        query_start (int): which query cds to start at
        query_domain_start (int): which query domain to start at
        target_index (dict): A dictionary containing tuples that correspond to the cds
        index of a domain, and the domain index of a domain
        target_start (int): which target cds to start at
        target_domain_start: which target domain to start at
        match (int): The score for a match between two domains
        mismatch (int): The penalty for a mismatch between two domains.
        gap (int): The penalty for a gap in the alignment
        max_match_dist (int): the maximum distance of a matching domain before it is
        considered a mismatch

    Returns:
        tuple: A tuple containing the query extension index on cds and domain level,
        target extension index on cds and domain level, and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0
    query_domain_exp = 0
    target_domain_exp = 0

    last_domain_idx = target_domain_start - 1

    for hsp_idx, hsp in enumerate(query_domains[query_domain_start:]):
        if hsp.domain not in target_index or not target_index[hsp.domain]:
            score += mismatch
            continue

        for dict_idx, target_idx in enumerate(target_index[hsp.domain]):
            cds_idx, domain_idx = target_idx

            match_after_lcs = True
            if domain_idx < (target_domain_start - 1):
                match_after_lcs = False
                continue

            # mismatch if the domain is too far away
            if abs(domain_idx - last_domain_idx) > max_match_dist:
                score += mismatch
                break

            # match. add score
            score += match

            # add gap penalties if the match is after the next position
            if domain_idx > last_domain_idx + 1:
                score += gap * (domain_idx - (last_domain_idx + 1))

            # subtract gap penalties if the match is before the current
            if domain_idx < last_domain_idx + 1:
                score -= gap

            if domain_idx > last_domain_idx:
                last_domain_idx = domain_idx

            # remove the current target index from the index
            target_index[hsp.domain].pop(dict_idx)

            break

        if not match_after_lcs:
            score += mismatch

        if score > max_score:
            max_score = score
            query_cds_idx = query_index[query_domain_start + hsp_idx]
            query_exp = query_cds_idx + 1 - query_start
            query_domain_exp = hsp_idx + 1
            if domain_idx + 1 - target_domain_start > target_domain_exp:
                target_exp = cds_idx + 1 - target_start
                target_domain_exp = domain_idx + 1 - target_domain_start

    return query_exp, query_domain_exp, target_exp, target_domain_exp, max_score


def score_extend_rev(
    query_domains: list[bs_hmm.HSP],
    query_index: dict[int, int],
    query_start: int,
    query_domain_start: int,
    target_index: dict[str, list[tuple[int, int]]],
    target_start: int,
    target_domain_start: int,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist: int,
) -> tuple[int, int, int, int, int]:
    """
    Calculate the score for extending a query sequence to a target sequence, in reverse

    TODO: This is a copy of score_extend. if possible, refactor to remove duplication

    Args:
        query_domains (list): A list of HSP objects representing the query record domains
        query_index (dict): A dictionary linking each domain index to its cds index
        query_start (int): which query cds to start at
        query_domain_start (int): which query domain to start at
        target_index (dict): A dictionary containing tuples that correspond to the cds
        index of a domain, and the domain index of a domain
        target_start (int): which target cds to start at
        target_domain_start: which target domain to start at
        match (int): The score for a match between two domains
        mismatch (int): The penalty for a mismatch between two domains.
        gap (int): The penalty for a gap in the alignment
        max_match_dist (int): the maximum distance of a matching domain before it is
        considered a mismatch

    Returns:
        tuple: A tuple containing the query extension index, target extension index,
        and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0
    query_domain_exp = 0
    target_domain_exp = 0

    # correct for inclusive start
    last_domain_idx = target_domain_start
    target_domain_start = target_domain_start - 1
    query_domain_start = query_domain_start - 1

    for hsp_idx, hsp in enumerate(query_domains[query_domain_start::-1]):
        if hsp.domain not in target_index or not target_index[hsp.domain]:
            score += mismatch
            continue

        for dict_idx, target_idx in enumerate(target_index[hsp.domain][::-1]):
            cds_idx, domain_idx = target_idx

            match_before_lcs = True
            if domain_idx > target_domain_start:
                match_before_lcs = False
                continue

            # mismatch if the domain is too far away
            if abs(domain_idx - last_domain_idx) > max_match_dist:
                score += mismatch
                break

            # match. add score
            score += match

            # add gap penalties if the match is after the next position
            if domain_idx < last_domain_idx - 1:
                score += gap * ((last_domain_idx - 1) - domain_idx)

            # subtract gap penalties if the match is before the current
            if domain_idx > last_domain_idx - 1:
                score -= gap

            if domain_idx < last_domain_idx:
                last_domain_idx = domain_idx

            # remove the current target index from the index
            target_index[hsp.domain].pop(dict_idx)

            break

        if not match_before_lcs:
            score += mismatch
            continue

        if score > max_score:
            max_score = score
            query_cds_idx = query_index[query_domain_start - hsp_idx]
            query_exp = query_start - query_cds_idx
            query_domain_exp = hsp_idx + 1
            if target_domain_start - domain_idx + 1 > target_domain_exp:
                target_exp = target_start - cds_idx
                target_domain_exp = target_domain_start - domain_idx + 1

    return query_exp, query_domain_exp, target_exp, target_domain_exp, max_score


def extend_greedy(pair: RecordPair) -> None:
    """Extends a comparable region in a greedy fashion

    This will extend the included set of cds in a pair as much as it can,
    based on the common domains found in the pair

    E.g. if we have the following two records:

    A: XAXXBXXXXCX
    B: XXXXXXXAXXXXXBCXXXXXXXXX

    The comparable region should be extended to the following:

    A: XAXXBXXXXCX
        [-------]
    B: XXXXXXXAXXXXXBCXXXXXXXXX
              [------]
    """
    logging.debug("before greedy extend:")
    logging.debug(pair.comparable_region)

    a_domains = list(pair.record_a.get_hsps())
    b_domains = list(pair.record_b.get_hsps())

    if pair.comparable_region.reverse:
        b_domains = b_domains[::-1]

    # index of domain to cds position
    a_index = get_target_indexes(a_domains)
    b_index = get_target_indexes(b_domains)

    common_domains = set(a_index.keys()).intersection(set(b_index.keys()))

    a_cds_min = len(list(pair.record_a.get_cds()))
    a_cds_max = 0
    b_cds_min = len(list(pair.record_b.get_cds()))
    b_cds_max = 0

    a_domain_min = len(a_domains)
    a_domain_max = 0
    b_domain_min = len(b_domains)
    b_domain_max = 0

    for domain in common_domains:
        for a_domain_index in a_index[domain]:
            a_cds_idx = a_domain_index[0]
            a_domain_idx = a_domain_index[1]

            a_cds_min = min(a_cds_min, a_cds_idx)
            a_cds_max = max(a_cds_max, a_cds_idx)
            a_domain_min = min(a_domain_min, a_domain_idx)
            a_domain_max = max(a_domain_max, a_domain_idx)

        for b_domain_index in b_index[domain]:
            b_cds_idx = b_domain_index[0]
            b_domain_idx = b_domain_index[1]

            b_cds_min = min(b_cds_min, b_cds_idx)
            b_cds_max = max(b_cds_max, b_cds_idx)
            b_domain_min = min(b_domain_min, b_domain_idx)
            b_domain_max = max(b_domain_max, b_domain_idx)

    pair.comparable_region.a_start = a_cds_min
    pair.comparable_region.a_stop = a_cds_max + 1
    pair.comparable_region.b_start = b_cds_min
    pair.comparable_region.b_stop = b_cds_max + 1

    pair.comparable_region.domain_a_start = a_domain_min
    pair.comparable_region.domain_a_stop = a_domain_max + 1
    pair.comparable_region.domain_b_start = b_domain_min
    pair.comparable_region.domain_b_stop = b_domain_max + 1

    logging.debug("after greedy extend:")
    logging.debug(pair.comparable_region)


def extend_simple_match(pair: RecordPair, match, gap):
    """Performs extension by first creating a simple match matrix, then
    performing a match/gap extension similar to legacy extension

    This method expects LCS to have been performed on the pair, and will
    do all four directions at once

    Args:
        pair: The record pair to extend
        match: The score for a match
        gap: The score for a gap
    """

    # I *think* the simplest way to do this is to create a list of domains with their source cds. for reasons.
    a_domains = []
    b_domains = []

    # so we'll do a loop through cds and through domains to keep track of everything
    for cds_idx, cds in enumerate(pair.record_a.get_cds_with_domains()):
        for domain in cds.hsps:
            a_domains.append((domain, cds_idx))

    for cds_idx, cds in enumerate(pair.record_b.get_cds_with_domains()):
        for domain in cds.hsps:
            b_domains.append((domain, cds_idx))

    # get the common domains
    common_domains = set([a[0] for a in a_domains]).intersection(
        set([b[0] for b in b_domains])
    )

    # TODO: this is obviously repetitive, but beyond extracting a function I'm
    # not sure how to make it better
    # for now I'd prefer this to be obvious and clear over looking good

    # a forward
    score = 0
    max_score = 0
    for i in range(pair.comparable_region.lcs_domain_a_stop, len(a_domains)):
        domain = a_domains[i][0]

        if domain not in common_domains:
            score += gap
        else:
            score += match

        if score > max_score:
            cds_idx = a_domains[i][1]
            max_score = score
            pair.comparable_region.a_stop = cds_idx + 1
            pair.comparable_region.domain_a_stop = i + 1

    # a reverse
    score = 0
    max_score = 0
    for i in range(pair.comparable_region.lcs_domain_a_start - 1, -1, -1):
        domain = a_domains[i][0]

        if domain not in common_domains:
            score += gap
        else:
            score += match

        if score > max_score:
            cds_idx = a_domains[i][1]
            max_score = score
            pair.comparable_region.a_start = cds_idx
            pair.comparable_region.domain_a_start = i

    # b forward
    score = 0
    max_score = 0
    for i in range(pair.comparable_region.lcs_domain_b_stop, len(b_domains)):
        domain = b_domains[i][0]

        if domain not in common_domains:
            score += gap
        else:
            score += match

        if score > max_score:
            cds_idx = b_domains[i][1]
            max_score = score
            pair.comparable_region.b_stop = cds_idx + 1
            pair.comparable_region.domain_b_stop = i + 1

    # b reverse
    score = 0
    max_score = 0
    for i in range(pair.comparable_region.lcs_domain_b_start - 1, -1, -1):
        domain = b_domains[i][0]

        if domain not in common_domains:
            score += gap
        else:
            score += match

        if score > max_score:
            cds_idx = b_domains[i][1]
            max_score = score
            pair.comparable_region.b_start = cds_idx
            pair.comparable_region.domain_b_start = i
