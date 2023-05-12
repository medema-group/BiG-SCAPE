"""Contains score extension code as it behaves in BiG-SCAPE 1.0"""

# from python
import logging

# from other modules
from src.parameters.constants import (
    EXPAND_MIN_LCS_LEN,
    EXPAND_GAP_SCORE,
    EXPAND_MATCH_SCORE,
    EXPAND_MISMATCH_SCORE,
)
from src.genbank import CDS
from src.comparison import ComparableRegion


def expand_glocal(
    comparable_region: ComparableRegion, min_match_len=EXPAND_MIN_LCS_LEN
) -> None:
    """Expand the comparable region on both sides using a simple scoring algorithm

    This assumes the initial a_len and b_len were set by LCS, and are the same.

    This is using a_len as LCS len
    """
    # exit early if the pair does not contain a core overlap of at least
    # EXPAND_MIN_LCS_LEN (3 in BiG-SCAPE 1.0) elements and no biosynthetic gene is found
    # in the LCS
    if comparable_region.a_stop < min_match_len:
        cds_stop = comparable_region.a_start + comparable_region.a_stop
        if not ComparableRegion.cds_range_contains_biosynthetic(
            comparable_region.pair.region_a, comparable_region.a_start, cds_stop, True
        ):
            logging.debug(
                (
                    "Skipping pair %s expansion - LCS len < %d and no core "
                    "biosynthetic genes found in LCS"
                ),
                comparable_region.pair,
                min_match_len,
            )
            return

    logging.debug("cr before expand: %s", str(comparable_region))

    expand_glocal_left(comparable_region)

    expand_glocal_right(comparable_region)

    logging.debug("cr after expand: %s", str(comparable_region))


def set_expansion_left(
    comparable_region: ComparableRegion, a_expansion: int, b_expansion: int
) -> None:
    """Set new start and length positions after expansion for the left side


    This function is similar to set_expansion_right, except that it only removes from
    the start of A

    This adds to the end of b only if B is reversed. otherwise it removes from
    the start position of B

    Illustrated:

    a_expansion: 3
    b_expansion: 3

    before:
    a: [-----LCS-----]
    b: [-----LCS-----]
    after, forward:
    a: [--<<<LCS-----]
    b: [--<<<LCS-----]
    after, reversed:
    a: [--<<<LCS-----]
    b: [-----LCS>>>--]

    Args:
        comparable_region (ComparableRegion): Comparable region for which to set
        expansion
        a_expansion (int): the length to expand A
        b_expansion (int): the length to expand B
    """
    comparable_region.a_start -= a_expansion
    comparable_region.b_start -= b_expansion


def expand_glocal_left(comparable_region: ComparableRegion) -> None:
    """Perform expansion on the left side of two regions where necessary and set new
    start and stop positions

    Args:
        comparable_region (ComparableRegion): comparable region for which to perform
        expansion
    """

    # this is where the fun begins
    # the legacy implementation of comparable region expansion first checks the number
    # of genes that are left of the current comparable region

    cds_list_a = comparable_region.pair.region_a.get_cds()

    a_left_stop = comparable_region.a_start - 1
    left_cds_a = cds_list_a[a_left_stop::-1]

    cds_list_b = comparable_region.pair.region_b.get_cds()

    if comparable_region.reverse:
        cds_list_b = cds_list_b[::-1]

    b_left_stop = comparable_region.b_start - 1
    left_cds_b = cds_list_b[b_left_stop::-1]

    # first check we do is to see which of the regions has more genes to the left

    # case 1: A has more genes to the left of LCS than B
    if len(left_cds_a) > len(left_cds_b):
        # in this case, we expand A based on B
        a_score, a_expansion = expand_score(left_cds_b, left_cds_a)
        # B is extended as far as it can be
        b_expansion = len(left_cds_b)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # case 2: B has more genes to the left of LCS than A
    if len(left_cds_b) > len(left_cds_a):
        # in this case, we expand B based on A
        b_score, b_expansion = expand_score(left_cds_a, left_cds_b)
        # A is extended as far as it can be
        a_expansion = len(left_cds_a)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # case 3: A and B have same number of genes left of LCS
    # first off, expand both
    a_score, a_expansion = expand_score(left_cds_b, left_cds_a)
    b_score, b_expansion = expand_score(left_cds_a, left_cds_b)

    # same score
    if a_score == b_score:
        # use A if it has a longer extension
        if a_expansion > b_expansion:
            set_expansion_left(comparable_region, a_expansion, a_expansion)
            return
        # otherwise just use B
        set_expansion_left(comparable_region, b_expansion, b_expansion)
        return

    # A has higher score
    if a_score > b_score:
        # ... use A
        set_expansion_left(comparable_region, a_expansion, a_expansion)
        return

    # only remaining case is B has higher score
    set_expansion_left(comparable_region, b_expansion, b_expansion)


def set_expansion_right(
    comparable_region: ComparableRegion, a_expansion: int, b_expansion: int
) -> None:
    """Set new start and length positions after expansion for the right side

    This function is similar to set_expansion_left, except that it only adds to the
    stop position of A.

    This adds to the stop position of b only if B was not reversed. otherwise it removes
    from the start position of B

    Illustrated:

    a_expansion: 3
    b_expansion: 3

    before:
    a: [-----LCS-----]
    b: [-----LCS-----]
    after, forward:
    a: [-----LCS>>>--]
    b: [-----LCS>>>--]
    after, reversed:
    a: [-----LCS>>>--]
    b: [--<<<LCS-----]


    Args:
        comparable_region (ComparableRegion): Comparable region for which to set
        expansion
        a_expansion (int): the length to expand A
        b_expansion (int): the length to expand B
    """
    comparable_region.a_stop += a_expansion
    comparable_region.b_stop += b_expansion


def expand_glocal_right(comparable_region: ComparableRegion) -> None:
    """Set new comparable region positions using glocal expand results

    Args:
    """

    # this is where the fun begins
    # the legacy implementation of comparable region expansion first checks the number
    # of genes that are right of the current comparable region

    cds_list_a = comparable_region.pair.region_a.get_cds()

    a_right_start = comparable_region.a_stop
    right_cds_a = cds_list_a[a_right_start:]

    cds_list_b = comparable_region.pair.region_b.get_cds()

    if comparable_region.reverse:
        cds_list_b = cds_list_b[::-1]

    b_right_start = comparable_region.b_stop
    right_cds_b = cds_list_b[b_right_start:]

    # first check we do is to see which of the regions has more genes to the left

    # case 1: A has more genes to the left of LCS than B
    if len(right_cds_a) > len(right_cds_b):
        # in this case, we expand A based on B
        a_score, a_expansion = expand_score(right_cds_a, right_cds_b)
        # B is extended as far as it can be
        b_expansion = len(right_cds_b)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # case 2: B has more genes to the left of LCS than A
    if len(right_cds_b) > len(right_cds_a):
        # in this case, we expand B based on A
        b_score, b_expansion = expand_score(right_cds_a, right_cds_b)
        # A is extended as far as it can be
        a_expansion = len(right_cds_a)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # case 3: A and B have same number of genes left of LCS
    # first off, expand both
    a_score, a_expansion = expand_score(right_cds_b, right_cds_a)
    b_score, b_expansion = expand_score(right_cds_a, right_cds_b)

    # same score
    if a_score == b_score:
        # use A if it has a longer extension
        if a_expansion > b_expansion:
            set_expansion_right(comparable_region, a_expansion, a_expansion)
            return
        # otherwise just use B
        set_expansion_right(comparable_region, b_expansion, b_expansion)
        return

    # A has higher score
    if a_score > b_score:
        # ... use A
        set_expansion_right(comparable_region, a_expansion, a_expansion)
        return

    # only remaining case is B has higher score
    set_expansion_right(comparable_region, b_expansion, b_expansion)


def expand_score(
    target: list[CDS],
    query: list[CDS],
    match_score: int = EXPAND_MATCH_SCORE,
    mismatch_score: int = EXPAND_MISMATCH_SCORE,
    gap_score: int = EXPAND_GAP_SCORE,
) -> tuple[int, int]:
    """TODO

    Args:
        target (lst[str]): the target to expand
        query (list[str]): the query upon which to base the expansion
        match_score (int): score to apply when a domain match is found
        mismatch_score (int): score to apply when a mismatch is found
        gap_score (int): score to apply when a gap is found

    Returns:
        tuple[int, int]: score, extension
    """

    max_score = 0
    extension = 0

    # TODO: exit loops early if there is no chance to recover to a good max score

    score = 0
    query_idx = 0
    skips = 0

    # expand right side.
    for target_idx, target_string in enumerate(target):
        # CDS with no domains
        if len(target_string.hsps) == 0:
            skips += 1

        # mismatch
        if target_string not in query[query_idx:]:
            score += mismatch_score
            continue

        # get the index of the match based upon where we left off last iteration
        match_pos = query[query_idx:].index(target_string)

        # because we are working from a slice of the query list, match_pos will tell us
        # how many elements we skip (gaps)
        score += match_score + match_pos * gap_score

        # update where we are at in the query
        query_idx += match_pos + 1

        # update max score
        if score >= max_score:
            max_score = score
            extension = target_idx + 1 + skips

    return max_score, extension