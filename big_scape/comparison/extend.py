"""Contains methods for extension of comparable regions with record pairs"""

# from python
import math
import logging

# from other modules
import big_scape.genbank as bs_genbank

# from this module
from .comparable_region import ComparableRegion


def reset(comparable_region: ComparableRegion) -> None:
    """Resets the expansion of a comparable region

    Args:
        comparable_region: The comparable region to reset
    """
    comparable_region.a_start = 0
    comparable_region.b_start = 0
    comparable_region.a_stop = len(
        comparable_region.pair.record_a.get_cds_with_domains()
    )
    comparable_region.b_stop = len(
        comparable_region.pair.record_b.get_cds_with_domains()
    )
    comparable_region.reverse = False


def check(
    comparable_region: ComparableRegion, min_len: int, biosynth_check: bool
) -> bool:
    """Checks if a comparable region should be reset after expansion

    returns true if either of the following conditions are met:
    - the comparable region contains a biosynthetic gene
    - the comparable region is longer than or equal to min_len

    Args:
        comparable_region: The comparable region to check
        min_len: The minimum length of the comparable region
        biosynth_check: Whether to check for biosynthetic genes within the comparable
            region
    """
    if biosynth_check:
        if ComparableRegion.cds_range_contains_biosynthetic(
            comparable_region.pair.record_a,
            comparable_region.a_start,
            comparable_region.a_stop,
            False,
            False,
        ):
            return True

        if ComparableRegion.cds_range_contains_biosynthetic(
            comparable_region.pair.record_b,
            comparable_region.b_start,
            comparable_region.b_stop,
            False,
            comparable_region.reverse,
        ):
            return True

    a_len = comparable_region.a_stop - comparable_region.a_start
    b_len = comparable_region.b_stop - comparable_region.b_start

    if a_len >= min_len and b_len >= min_len:
        return True

    return False


def extend(
    comparable_region: ComparableRegion,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist_perc: float,
) -> None:
    """Expands a comparable region

    This will expand the included set of cds in a comparable region based on a scoring
    mechanism. If the pair in the comparable region consists of protoclusters, the
    this will not be limited to the bounds of those protoclusters

    NOTE: The start-stop slices returned from this function correspond to CDS slices,
    without domains only.

    Args:
        comparable_region: The comparable region to expand
        match: The score for a match
        mismatch: The score for a mismatch
        gap: The score for a gap
        max_match_dist_perc: The maximum distance of a matching domain before it is
    """

    logging.debug("before extend:")
    logging.debug(comparable_region)

    # get the cds lists
    a_cds = comparable_region.pair.record_a.get_cds_with_domains(True)
    b_cds = comparable_region.pair.record_b.get_cds_with_domains(True)

    # TODO: base on domains not cds
    a_max_dist = math.floor(len(a_cds) * max_match_dist_perc)
    b_max_dist = math.floor(len(b_cds) * max_match_dist_perc)

    # reverse b if necessary. This might be true after LCS
    if comparable_region.reverse:
        b_cds = b_cds[::-1]

    # we will try the following approach:
    # the shorter cds is the query
    # the longer cds is the target
    # generate a dictionary of domain positions in the target
    # we try to expand the target by finding matching domains in the query
    # when we do not find a match, this counts as a mismatch and we add a penalty
    # when we find a match, add gap penalties for the number of gaps we have to insert
    # if we find a match before the current position, subtract a gap penalty
    # we will only match when the position from the current extension is below a cutoff
    # length

    if len(a_cds) > len(b_cds):
        query = b_cds
        query_start = comparable_region.b_start
        query_stop = comparable_region.b_stop
        target = a_cds
        target_start = comparable_region.a_start
        target_stop = comparable_region.a_stop
        max_match_dist = a_max_dist
    else:
        query = a_cds
        query_start = comparable_region.a_start
        query_stop = comparable_region.a_stop
        target = b_cds
        target_start = comparable_region.b_start
        target_stop = comparable_region.b_stop
        max_match_dist = b_max_dist

    if target_stop != len(target) and query_stop != len(query):
        # generate an index of domain positions in the target
        # the lists in this index will be sorted by position, in ascending order
        target_index, target_cds_to_domain_index = get_target_indexes(target)

        query_exp, target_exp, score = score_extend(
            query,
            query_stop,
            target_index,
            target_stop,
            target_cds_to_domain_index[target_stop],
            match,
            mismatch,
            gap,
            max_match_dist,
        )

        # set the new start and stop positions
        if len(a_cds) > len(b_cds):
            comparable_region.b_stop += query_exp
            comparable_region.a_stop += target_exp
        else:
            comparable_region.a_stop += query_exp
            comparable_region.b_stop += target_exp

    if target_start != 0 and query_start != 0:
        # reset indexes as they could have been modified in forward extend
        target_index, target_cds_to_domain_index = get_target_indexes(target)

        query_exp, target_exp, score = score_extend_rev(
            query,
            query_start,
            target_index,
            target_start,
            target_cds_to_domain_index[target_start],
            match,
            mismatch,
            gap,
            max_match_dist,
        )

        # expand left
        if len(a_cds) > len(b_cds):
            comparable_region.b_start -= query_exp
            comparable_region.a_start -= target_exp
        else:
            comparable_region.a_start -= query_exp
            comparable_region.b_start -= target_exp

    logging.debug("after extend:")
    logging.debug(comparable_region)


def get_target_indexes(
    target: list[bs_genbank.CDS],
) -> tuple[dict[str, list[tuple[int, int]]], dict[int, int]]:
    """Generate two indexes of domain cds positions in the target

    First dictionary that is returned will have the following structure:
    {
        domain: [(cds_idx, domain_idx), ...]
    }

    All indexes are sorted by position, in ascending order

    Second dictionary stores the first domain idx for each cds: {cds_idx: domain_idx}

    Args:
        target (list[bs_genbank.CDS]): The target cds list

    Returns:
        tuple[dict[str, list[tuple[int, int]]], dict[int, int]]: The target index
        dictionaries
    """

    target_index: dict[str, list[tuple[int, int]]] = {}
    target_cds_to_domain_index: dict[int, int] = {}
    domain_idx = 0
    for cds_idx, cds in enumerate(target):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                target_index[hsp.domain] = []

            target_index[hsp.domain].append((cds_idx, domain_idx))

            if cds_idx not in target_cds_to_domain_index:
                target_cds_to_domain_index[cds_idx] = domain_idx

            domain_idx += 1
    return target_index, target_cds_to_domain_index


def score_extend(
    query: list[bs_genbank.CDS],
    query_start: int,
    target_index: dict[str, list[tuple[int, int]]],
    target_start: int,
    target_domain_start: int,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist: int,
) -> tuple[int, int, int]:
    """
    Calculate the score for extending a query sequence to a target sequence.

    TODO: This is a copy of score_extend_rev. if possible, refactor to remove
    duplication

    Args:
        query (list): A list of Cds objects representing the full query record cds
        query_start (int): which query cds to start at
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
        tuple: A tuple containing the query expansion index, target expansion index,
        and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0

    last_domain_idx = target_domain_start - 1

    for query_idx, cds in enumerate(query[query_start:]):
        for hsp in cds.hsps:
            if hsp.domain not in target_index or not target_index[hsp.domain]:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain]):
                cds_idx, domain_idx = target_idx

                match_after_lcs = True
                if cds_idx < target_start:
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
                query_exp = query_idx + 1
                if cds_idx + 1 - target_start > target_exp:
                    target_exp = cds_idx + 1 - target_start

    return query_exp, target_exp, max_score


def score_extend_rev(
    query: list[bs_genbank.CDS],
    query_start: int,
    target_index: dict[str, list[tuple[int, int]]],
    target_start: int,
    target_domain_start: int,
    match: int,
    mismatch: int,
    gap: int,
    max_match_dist: int,
) -> tuple[int, int, int]:
    """
    Calculate the score for extending a query sequence to a target sequence, in reverse

    TODO: This is a copy of score_extend. if possible, refactor to remove duplication

    Args:
        query (list): A list of Cds objects representing the full query record cds
        query_start (int): which query cds to start at
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
        tuple: A tuple containing the query expansion index, target expansion index,
        and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0

    # correct for inclusive start
    last_domain_idx = target_domain_start
    target_start = target_start - 1
    query_start = query_start - 1

    for query_idx, cds in enumerate(query[query_start::-1]):
        for hsp in cds.hsps[::-1]:
            if hsp.domain not in target_index or not target_index[hsp.domain]:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain][::-1]):
                cds_idx, domain_idx = target_idx

                match_before_lcs = True
                if cds_idx > target_start:
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
                query_exp = query_idx + 1
                if target_start - cds_idx + 1 > target_exp:
                    target_exp = target_start - cds_idx + 1

    return query_exp, target_exp, max_score
