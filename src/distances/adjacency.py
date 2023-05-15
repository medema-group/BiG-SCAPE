"""Contains code to calculate adjacency indexes of set pairs"""


# from other modules
from src.comparison import BGCPair


def calc_ai_lists(list_a: list, list_b: list):
    """Calculates the adjacency index of two lists, which is the Jaccard index of sets
    of neighbouring items in two sorted lists.

    e.g.

    A: [a,b,c]
    B: [b,c,e]

    a neighbours: [ab, bc]
    b neighbours: [bc, ce]

    intersect: [bc] (1)
    union: [ab, bc, ce] (3)
    jaccard: 1/3
    """
    if len(list_a) < 2:
        return 0.0

    if len(list_b) < 2:
        return 0.0

    set_a = set()
    set_b = set()

    for idx_a in range(len(list_a) - 1):
        start = idx_a
        stop = start + 2
        neighbours_a = tuple(sorted(list_a[start:stop]))

        set_a.add(neighbours_a)

    for idx_b in range(len(list_b) - 1):
        start = idx_b
        stop = start + 2
        neighbours_b = tuple(sorted(list_b[start:stop]))

        set_b.add(neighbours_b)

    intersect = set_a & set_b
    union = set_a | set_b

    adjacency_index = len(intersect) / len(union)

    return adjacency_index


def calc_ai_pair(bgc_pair: BGCPair) -> float:
    """Calculate the adjacency index of a pair of BGCs

    Args:
        pair (BGCPair): Pair of BGCs for which to calculate the AI

    Returns:
        float: adjacency index
    """
    a_domains, b_domains = bgc_pair.comparable_region.get_domain_lists()

    a_domain_names = [hsp.domain for hsp in a_domains]
    b_domain_names = [hsp.domain for hsp in b_domains]

    return calc_ai_lists(a_domain_names, b_domain_names)
