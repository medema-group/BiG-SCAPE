"""Contains tests for legacy methods that intend to mimic BiG-SCAPE 1.0 behaviors"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import CDS
from src.hmm import HSP
from src.comparison.legacy_extend import expand_score


class TestLegacyComparableRegion(TestCase):
    """Contains tests for legacy comparable region methods"""

    def test_legacy_expand(self):
        """Tests whether expand_score correctly returns a length and score for a given
        sequence of strings
        """
        # LCS is CDE, valid expansion is ACDE
        target_domains = ["A", "C", "D", "E"]
        query_domains = ["A", "B", "C", "D", "E"]

        expand_start = 3

        target_cds_list = []
        for target_domain in target_domains:
            target_cds = CDS(0, 10)
            target_cds.hsps.add(HSP(target_cds, target_domain, 100, 0, 10))
            target_cds_list.append(target_cds)

        query_cds_list = []
        for query_domain in query_domains:
            query_cds = CDS(0, 10)
            query_cds.hsps.add(HSP(query_cds, query_domain, 100, 0, 10))
            query_cds_list.append(query_cds)

        # gap = 1 * -2 = -2
        # match = 1 * 5 = 5
        # score = 5 - 2 = 3
        # expansion = 1
        expected_results = (3, 1)

        actual_results = expand_score(
            target_cds_list[expand_start:], query_cds_list[expand_start:], 1, 2, 3
        )

        self.assertEqual(expected_results, actual_results)
