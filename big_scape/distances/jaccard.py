"""Contains code to calculate the Jaccard index of two BGCs"""

# from python
from typing import Any

# from other modules
from big_scape.comparison import RecordPair


def calc_jaccard_sets(set_a: set[Any], set_b: set[Any]) -> float:
    """Calculates the Jaccard index between two sets of objects

    Args:
        set_a (set[Any]): First set
        set_b (set[Any]): Second set

    Returns:
        float: The Jaccard index of these two sets
    """
    if len(set_a) == 0 or len(set_b) == 0:
        return 0

    return len(set_a & set_b) / len(set_a | set_b)


def calc_jaccard_pair(bgc_pair: RecordPair, cache=True) -> float:
    """Generates the Jaccard index between a pair of BGCs

    Args:
        bgc_pair (BGCPair): The pair of BGCs to calculate the Jaccard index for

    Returns:
        float: Jaccard index
    """
    a_cds, b_cds = bgc_pair.get_domain_sets(cache=cache)
    return calc_jaccard_sets(a_cds, b_cds)
