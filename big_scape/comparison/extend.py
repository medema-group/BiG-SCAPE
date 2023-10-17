"""Contains methods for extension of comparable regions with record pairs"""

# from python


# from this modile
from .comparable_region import ComparableRegion


def reset_expansion(comparable_region: ComparableRegion) -> None:
    """Resets the expansion of a comparable region

    Args:
        comparable_region: The comparable region to reset
    """
    comparable_region.a_start = 0
    comparable_region.b_start = 0
    comparable_region.a_stop = len(comparable_region.pair.region_a.get_cds())
    comparable_region.b_stop = len(comparable_region.pair.region_b.get_cds())
    comparable_region.reverse = False


def check_expand(comparable_region: ComparableRegion, min_len=3) -> bool:
    """Checks if a comparable region should be reset after expansion

    returns true if either of the following conditions are met:
    - the comparable region contains a biosynthetic gene
    - the comparable region is longer than or equal to min_len

    Args:
        comparable_region: The comparable region to check
    """
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


def check_lcs(comparable_region: ComparableRegion) -> None:
    """Checks if a comparable region should be expanded

    Args:
        comparable_region: The comparable region to check
    """
    pass


def expand(comparable_region: ComparableRegion) -> None:
    """Expands a comparable region

    Args:
        comparable_region: The comparable region to expand
    """
    pass
