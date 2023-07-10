"""Contains score extension code as it behaves in BiG-SCAPE 1.0"""

# from python
import logging

# from other modules
from src.parameters.constants import (
    MIN_LCS_LEN,
    MIN_EXPAND_LEN,
    EXPAND_GAP_SCORE,
    EXPAND_MATCH_SCORE,
    EXPAND_MISMATCH_SCORE,
)
from src.genbank import CDS
from src.comparison import ComparableRegion, BGCPair


def reset_expansion(
    comparable_region: ComparableRegion, a_start=0, a_stop=None, b_start=0, b_stop=None
) -> None:
    """Resets the comparable region starts and stops

    If no arguments beyond comparable region are given, resets to the full range
    """
    a_gbk = comparable_region.pair.region_a.parent_gbk
    b_gbk = comparable_region.pair.region_b.parent_gbk

    if a_gbk is None or b_gbk is None:
        return

    if a_stop is None:
        a_stop = len(a_gbk.genes) + 1

    if b_stop is None:
        b_stop = len(b_gbk.genes) + 1

    comparable_region.a_start = a_start
    comparable_region.b_start = b_start

    comparable_region.a_stop = a_stop
    comparable_region.b_stop = b_stop

    comparable_region.domain_lists = None
    comparable_region.domain_sets = None


def expand_glocal(
    comparable_region: ComparableRegion,
    min_lcs_len=MIN_LCS_LEN,
    min_expand_len=MIN_EXPAND_LEN,
) -> None:
    """Expand the comparable region on both sides using a simple scoring algorithm

    This assumes the initial a_len and b_len were set by LCS, and are the same.

    This is using a_len as LCS len
    """
    # exit early if the pair does not contain a core overlap of at least
    # EXPAND_MIN_LCS_LEN (3 in BiG-SCAPE 1.0) elements and no biosynthetic gene is found
    # in the LCS
    # TODO: are all relevant checks here?
    match_len = comparable_region.a_stop - comparable_region.a_start
    if match_len < min_lcs_len:
        if not ComparableRegion.cds_range_contains_biosynthetic(
            comparable_region.pair.region_a,
            comparable_region.a_start,
            comparable_region.a_stop,
            True,
        ):
            logging.debug(
                (
                    "Skipping pair %s left expansion - LCS len < %d and no core "
                    "biosynthetic genes found in LCS"
                ),
                comparable_region.pair,
                min_lcs_len,
            )
            reset_expansion(comparable_region)
            return

    logging.debug("cr before expand: %s", str(comparable_region))

    expand_glocal_left(comparable_region)

    expand_glocal_right(comparable_region)

    logging.debug("cr after expand: %s", str(comparable_region))


def check_expand(
    comparable_region: ComparableRegion, min_expand_len=MIN_EXPAND_LEN
) -> bool:
    """Returns True if the expansion is valid. Returns False if the expansion should be reset"""

    # final checks: did we expand enough?
    a_start = comparable_region.a_start
    a_stop = comparable_region.a_stop
    cds_list_a = comparable_region.pair.region_a.get_cds_with_domains()[a_start:a_stop]
    expansion_len_a = len([cds for cds in cds_list_a if len(cds.hsps) > 0])

    b_start = comparable_region.b_start
    b_stop = comparable_region.b_stop
    cds_list_b = comparable_region.pair.region_b.get_cds_with_domains()[b_start:b_stop]
    expansion_len_b = len([cds for cds in cds_list_b if len(cds.hsps) > 0])

    if min(expansion_len_a, expansion_len_b) < min_expand_len:
        return False

    # do both slices contain a biosynthetic gene?
    if not ComparableRegion.cds_range_contains_biosynthetic(
        comparable_region.pair.region_a,
        comparable_region.a_start,
        comparable_region.a_stop,
        end_inclusive=True,
    ):
        return False

    if not ComparableRegion.cds_range_contains_biosynthetic(
        comparable_region.pair.region_b,
        comparable_region.b_start,
        comparable_region.b_stop,
        end_inclusive=True,
        reverse=comparable_region.reverse,
    ):
        return False

    return True


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

    cds_list_a = comparable_region.pair.region_a.get_cds_with_domains()

    a_left_stop = comparable_region.a_start
    left_cds_a = list(reversed(cds_list_a[:a_left_stop]))
    left_domain_cds_a = len([cds for cds in left_cds_a if len(cds.hsps) > 0])

    cds_list_b = comparable_region.pair.region_b.get_cds_with_domains(
        reverse=comparable_region.reverse
    )

    b_left_stop = comparable_region.b_start
    left_cds_b = list(reversed(cds_list_b[:b_left_stop]))
    left_domain_cds_b = len([cds for cds in left_cds_b if len(cds.hsps) > 0])

    # first check we do is to see which of the regions has more genes to the left

    # case 1: A has more genes to the left of LCS than B
    if left_domain_cds_a > left_domain_cds_b:
        # in this case, we expand A based on B
        a_score, a_expansion = expand_score(left_cds_a, left_cds_b)
        # B is extended as far as it can be
        b_expansion = len(left_cds_b)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # case 2: B has more genes to the left of LCS than A
    if left_domain_cds_b > left_domain_cds_a:
        # in this case, we expand B based on A
        b_score, b_expansion = expand_score(left_cds_b, left_cds_a)
        # A is extended as far as it can be
        a_expansion = len(left_cds_a)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # case 3: A and B have same number of genes left of LCS
    # first off, expand both
    a_score, a_expansion = expand_score(left_cds_a, left_cds_b)
    b_score, b_expansion = expand_score(left_cds_b, left_cds_a)

    # same score
    if a_score == b_score:
        # use A if it has a longer extension
        if a_expansion > b_expansion:
            b_expansion = len(left_cds_b)
            set_expansion_left(comparable_region, a_expansion, b_expansion)
            return
        # otherwise just use B
        a_expansion = len(left_cds_a)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # A has higher score
    if a_score > b_score:
        # ... use A
        b_expansion = len(left_cds_b)
        set_expansion_left(comparable_region, a_expansion, b_expansion)
        return

    # only remaining case is B has higher score
    a_expansion = len(left_cds_a)
    set_expansion_left(comparable_region, a_expansion, b_expansion)


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

    cds_list_a = comparable_region.pair.region_a.get_cds_with_domains()

    a_right_start = comparable_region.a_stop
    right_cds_a = cds_list_a[a_right_start:]
    right_domain_cds_a = len([cds for cds in right_cds_a if len(cds.hsps) > 0])

    cds_list_b = comparable_region.pair.region_b.get_cds_with_domains(
        reverse=comparable_region.reverse
    )

    b_right_start = comparable_region.b_stop
    right_cds_b = cds_list_b[b_right_start:]
    right_domain_cds_b = len([cds for cds in right_cds_b if len(cds.hsps) > 0])

    # first check we do is to see which of the regions has more genes to the left

    # case 1: A has more genes to the left of LCS than B
    if right_domain_cds_a > right_domain_cds_b:
        # in this case, we expand A based on B
        a_score, a_expansion = expand_score(right_cds_a, right_cds_b)
        # B is extended as far as it can be
        b_expansion = len(right_cds_b)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # case 2: B has more genes to the left of LCS than A
    if right_domain_cds_b > right_domain_cds_a:
        # in this case, we expand B based on A
        b_score, b_expansion = expand_score(right_cds_b, right_cds_a)
        # A is extended as far as it can be
        a_expansion = len(right_cds_a)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # case 3: A and B have same number of genes left of LCS
    # first off, expand both
    a_score, a_expansion = expand_score(right_cds_a, right_cds_b)
    b_score, b_expansion = expand_score(right_cds_b, right_cds_a)

    # same score
    if a_score == b_score:
        # use A if it has a longer extension
        if a_expansion > b_expansion:
            b_expansion = len(right_cds_b)
            set_expansion_right(comparable_region, a_expansion, b_expansion)
            return
        # otherwise just use B
        a_expansion = len(right_cds_a)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # A has higher score
    if a_score > b_score:
        # ... use A
        b_expansion = len(right_cds_b)
        set_expansion_right(comparable_region, a_expansion, b_expansion)
        return

    # only remaining case is B has higher score
    a_expansion = len(right_cds_a)
    set_expansion_right(comparable_region, a_expansion, b_expansion)


def expand_score(
    query: list[CDS],
    target: list[CDS],
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
    score = 0
    extension = 0

    # TODO: exit loops early if there is no chance to recover to a good max score

    target_idx = 0

    # expand right side.
    for query_idx, query_string in enumerate(query):
        # mismatch
        if query_string not in target[target_idx:]:
            score += mismatch_score
            continue

        # get the index of the match based upon where we left off last iteration
        match_pos = target[target_idx:].index(query_string)

        # because we are working from a slice of the query list, match_pos will tell us
        # how many elements we skip (gaps)
        score += match_score + match_pos * gap_score

        # update where we are at in the query
        target_idx += match_pos + 1

        # update max score
        if score >= max_score:
            max_score = score
            extension = query_idx + 1

    return max_score, extension


def legacy_needs_extend(
    pair: BGCPair, alignment_mode: str, extend_slice_cutoff=MIN_LCS_LEN
):
    """Returns False if:

    - alignment_mode is global
    - alignment mode is auto, and:
        pair.region_a and pair_region_b are both not on a contig edge
    - region_a does not contain a biosynthetic gene and:
        comparable region length of a is less than extend_slice_cutoff
        (at this point, this should be the LCS length)

    """

    if alignment_mode == "global":
        return False

    if alignment_mode == "auto" and not (
        pair.region_a.contig_edge and pair.region_b.contig_edge
    ):
        return False

    lcs_extend_len = pair.comparable_region.a_stop - 1 - pair.comparable_region.a_start
    if (
        lcs_extend_len < extend_slice_cutoff
        and not ComparableRegion.cds_range_contains_biosynthetic(
            pair.region_a,
            pair.comparable_region.a_start,
            pair.comparable_region.a_stop,
            end_inclusive=True,  # technically wrong, but 1.0 behavior
            reverse=False,
        )
    ):
        return False

    return True
