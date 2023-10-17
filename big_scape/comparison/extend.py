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


def check_expand(comparable_region: ComparableRegion) -> None:
    """Checks if a comparable region should be reset after expansion, e.g. if lengths
    are too short

    Args:
        comparable_region: The comparable region to check
    """
    pass


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
