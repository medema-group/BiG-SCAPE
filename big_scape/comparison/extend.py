"""Contains methods for extension of comparable regions with record pairs"""

# from python
import math

# from other modules
import big_scape.genbank as bs_genbank

# from this module
import logging
from .comparable_region import ComparableRegion


def reset(comparable_region: ComparableRegion) -> None:
    """Resets the expansion of a comparable region

    Args:
        comparable_region: The comparable region to reset
    """
    comparable_region.a_start = 0
    comparable_region.b_start = 0
    comparable_region.a_stop = len(
        comparable_region.pair.region_a.get_cds_with_domains()
    )
    comparable_region.b_stop = len(
        comparable_region.pair.region_b.get_cds_with_domains()
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
            comparable_region.pair.region_a,
            comparable_region.a_start,
            comparable_region.a_stop,
            False,
            False,
        ):
            return True

        if ComparableRegion.cds_range_contains_biosynthetic(
            comparable_region.pair.region_b,
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
    a_cds = comparable_region.pair.region_a.get_cds_with_domains(True)
    b_cds = comparable_region.pair.region_b.get_cds_with_domains(True)

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
        max_match_dist = a_max_dist
    else:
        query = a_cds
        query_start = comparable_region.a_start
        query_stop = comparable_region.a_stop
        target = b_cds
        max_match_dist = b_max_dist

    # generate an index of domain positions in the target
    # the lists in this index will be sorted by position, in ascending order
    target_index = get_target_indexes(target)

    query_exp, target_exp, score = score_extend(
        query, query_stop, target_index, match, mismatch, gap, max_match_dist
    )

    # set the new start and stop positions
    if len(a_cds) > len(b_cds):
        comparable_region.b_stop += query_exp
        comparable_region.a_stop += target_exp
    else:
        comparable_region.a_stop += query_exp
        comparable_region.b_stop += target_exp

    query_exp, target_exp, score = score_extend_rev(
        query, query_start, target_index, match, mismatch, gap, max_match_dist
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
) -> dict[str, list[tuple[int, int]]]:
    """Generate an index of domain cds positions in the target

    The dictionary that is returned will have the following structure:
    {
        domain: [(cds_idx, domain_idx), ...]
    }

    All indexes are sorted by position, in ascending order

    Args:
        target (list[bs_genbank.CDS]): The target cds list

    Returns:
        dict[str, list[tuple[int, int]]]: The target index dictionary
    """

    target_index: dict[str, list[tuple[int, int]]] = {}
    domain_idx = 0
    for cds_idx, cds in enumerate(target):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                target_index[hsp.domain] = []

            target_index[hsp.domain].append((cds_idx, domain_idx))
            domain_idx += 1

    return target_index


def score_extend(
    query, query_start, target_index, match, mismatch, gap, max_match_dist
):
    """
    Calculate the score for extending a query sequence to a target sequence.

    TODO: max_match_dist
    TODO: This is a copy of score_extend_rev. if possible, refactor to remove
    duplication

    Args:
        query (list): A list of Cds objects representing the full query record cds
        query_start (int): which cds to start at
        target_index (dict): A dictionary containing tuples that correspond to the cds
        index of a domain, and the domain index of a domain
        match (int): The score for a match between two domains
        mismatch (int): The penalty for a mismatch between two domains.
        gap (int): The penalty for a gap in the alignment
        max_match_dist (int): the maximum distance of a matching domain before it is
        considered a mismatch

    Returns:
        tuple: A tuple containing the target expansion index, query expansion index, and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0

    for query_idx, cds in enumerate(query[query_start:]):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain]):
                cds_idx, domain_idx = target_idx

                if domain_idx < query_start:
                    continue

                # match. add score
                score += match

                # add gap penalties if the match is after the next position
                if domain_idx > target_exp:
                    score += gap * (domain_idx - (target_exp))

                # subtract gap penalties if the match is before the current
                if domain_idx < target_exp:
                    score -= gap

                # remove the current target index from the index
                target_index[hsp.domain].pop(dict_idx)

                break

            if score > max_score:
                max_score = score
                query_exp = query_idx + 1
                if cds_idx + 1 > target_exp:
                    target_exp = cds_idx + 1

    return target_exp, query_exp, max_score


def score_extend_rev(
    query, query_start, target_index, match, mismatch, gap, max_match_dist
):
    """
    Calculate the score for extending a query sequence to a target sequence, in reverse

    TODO: This is a copy of score_extend. if possible, refactor to remove duplication

    Args:
        query (list): A list of Cds objects representing the full query record cds
        query_start (int): which cds to start at
        target_index (dict): A dictionary containing tuples that correspond to the cds
        index of a domain, and the domain index of a domain
        match (int): The score for a match between two domains
        mismatch (int): The penalty for a mismatch between two domains.
        gap (int): The penalty for a gap in the alignment
        max_match_dist (int): the maximum distance of a matching domain before it is
        considered a mismatch

    Returns:
        tuple: A tuple containing the target expansion index, query expansion index,
        and the maximum score.
    """
    score = 0
    max_score = 0
    target_exp = 0
    query_exp = 0

    for query_idx, cds in enumerate(query[query_start::-1]):
        query_idx = len(query) - query_idx - 1
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain][::-1]):
                cds_idx, domain_idx = target_idx

                if domain_idx < query_start:
                    continue

                # match. add score
                score += match

                # add gap penalties if the match is after the next position
                if domain_idx > target_exp:
                    score += gap * (domain_idx - (target_exp))

                # subtract gap penalties if the match is before the current
                if domain_idx < target_exp:
                    score -= gap

                # remove the current target index from the index
                target_index[hsp.domain].pop(dict_idx)

                break

            if score > max_score:
                max_score = score
                query_exp = query_idx + 1
                if cds_idx + 1 > target_exp:
                    target_exp = cds_idx + 1

    return target_exp, query_exp, max_score
