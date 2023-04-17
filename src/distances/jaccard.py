"""Contains code to calculate the Jaccard index of two BGCs"""


# from other modules
from typing import Any
from src.comparison import BGCPair


def calculate_jaccard_sets(set_a: set[Any], set_b: set[Any]) -> float:
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


def calculate_jaccard_pair(bgc_pair: BGCPair) -> float:
    """Generates the Jaccard index between a pair of BGCs

    Args:
        bgc_pair (BGCPair): The pair of BGCs to calculate the Jaccard index for

    Returns:
        float: Jaccard index
    """
    a_genes = bgc_pair.region_a.get_cds()
    b_genes = bgc_pair.region_b.get_cds()

    a_domains: set[str] = set()
    b_domains: set[str] = set()

    for a_gene in a_genes:
        for a_hsp in a_gene.hsps:
            a_domains.add(a_hsp.domain)

    for b_gene in b_genes:
        for b_hsp in b_gene.hsps:
            b_domains.add(b_hsp.domain)

    return calculate_jaccard_sets(a_domains, b_domains)
