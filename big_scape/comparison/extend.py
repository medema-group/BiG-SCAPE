"""Contains methods for extension of comparable regions with record pairs"""

# from python


# from this modile
import logging
from .comparable_region import ComparableRegion


def reset(comparable_region: ComparableRegion) -> None:
    """Resets the expansion of a comparable region

    Args:
        comparable_region: The comparable region to reset
    """
    comparable_region.a_start = 0
    comparable_region.b_start = 0
    comparable_region.a_stop = len(comparable_region.pair.region_a.get_cds())
    comparable_region.b_stop = len(comparable_region.pair.region_b.get_cds())
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


def extend(comparable_region: ComparableRegion) -> None:
    """Expands a comparable region

    This will expand the included set of cds in a comparable region based on a scoring
    mechanism. If the pair in the comparable region consists of protoclusters, the
    this will not be limited to the bounds of those protoclusters

    Args:
        comparable_region: The comparable region to expand
    """
    logging.info(comparable_region.pair)
    logging.info("before:")
    logging.info(comparable_region)

    # get the cds lists
    a_cds = comparable_region.pair.region_a.get_cds(True)
    cds_b = comparable_region.pair.region_b.get_cds(True)

    # reverse b if necessary. This might be true after LCS
    if comparable_region.reverse:
        cds_b = cds_b[::-1]

    # we will try the following approach:
    # the shorter cds is the query
    # the longer cds is the target
    # generate a dictionary of domain positions in the target
    # we try to expand the target by finding matching domains in the query
    # when we do not find a match, this counts as a mismatch and we add a penalty
    # when we find a match, add gap penalties for the number of gaps we have to insert
    # if we find a match before the current position, subtract a gap penalty
    # TODO: how do we handle cases where we can have a match in other end of the cds

    if len(a_cds) > len(cds_b):
        query = cds_b
        query_start = comparable_region.b_start
        query_stop = comparable_region.b_stop
        target = a_cds
    else:
        query = a_cds
        query_start = comparable_region.a_start
        query_stop = comparable_region.a_stop
        target = cds_b

    # generate an index of domain positions in the target
    # the lists in this index will be sorted by position
    target_index: dict[str, list[int]] = {}
    for i, cds in enumerate(target):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                target_index[hsp.domain] = []
            target_index[hsp.domain].append(i)

    match = 5
    mismatch = -3
    gap = -2

    # start on the right side
    max_score = 0
    score = 0
    target_exp = 0
    query_exp = 0

    for query_idx, cds in enumerate(query[query_stop:]):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain]):
                if target_idx < query_stop:
                    continue

                # match. add score
                score += match

                # add gap penalties if the match is after the current position
                if target_idx > target_exp:
                    score += gap * (target_idx - target_exp)

                # subtract gap penalties if the match is before the current position
                if target_idx < target_exp:
                    score -= gap

                # remove the current target index from the index
                target_index[hsp.domain].pop(dict_idx)

                break

            if score > max_score:
                max_score = score
                query_exp = query_idx
                target_exp = target_idx

    # set the new start and stop positions
    if len(a_cds) > len(cds_b):
        comparable_region.b_stop += query_exp
        comparable_region.a_stop += target_exp
    else:
        comparable_region.a_stop += query_exp
        comparable_region.b_stop += target_exp

    # left side
    score = 0
    target_exp = 0
    query_exp = 0

    # flip the sequences
    query = query[::-1]
    target = target[::-1]

    # flip the start
    query_start = len(query) - query_start

    for query_idx, cds in enumerate(query[query_start:]):
        for hsp in cds.hsps:
            if hsp.domain not in target_index:
                score += mismatch
                continue

            for dict_idx, target_idx in enumerate(target_index[hsp.domain]):
                if target_idx > query_start:
                    continue

                # match. add score
                score += match

                # add gap penalties if the match is after the current position
                if target_idx > target_exp:
                    score += gap * (target_idx - target_exp)

                # subtract gap penalties if the match is before the current position
                if target_idx < target_exp:
                    score -= gap

                # remove the current target index from the index
                target_index[hsp.domain].pop(dict_idx)

                break

            if score > max_score:
                max_score = score
                query_exp = query_idx
                target_exp = target_idx

    # expand left
    if len(a_cds) > len(cds_b):
        comparable_region.b_start -= query_exp
        comparable_region.a_start -= target_exp
    else:
        comparable_region.a_start -= query_exp
        comparable_region.b_start -= target_exp

    logging.info("after:")
    logging.info(comparable_region)
    logging.info("")
