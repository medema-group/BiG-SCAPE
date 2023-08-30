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
        # LCS is CDE, valid expansion is CDEGH
        query_expand_start = 5
        query_domains = [
            "A",
            "B",
            "C",  # LCS
            "D",  # LCS
            "E",  # LCS
            "F",
            "H",
            "I",
            "J",
        ]

        target_expand_start = 4
        target_domains = [
            "A",  # before LCS
            # gap, before LCS
            "C",  # LCS
            "D",  # LCS
            "E",  # LCS
            # gap -2
            "G",  # mismatch -3
            "H",  # match +5
            # gap -2
            "J",  # match +5
        ]

        # 0 - 3 + 5 - 2 + 5 - 2 = 3
        # expands 3 positions
        expected_results = (3, 3)

        target_cds_list = []
        for target_domain in target_domains:
            target_cds = CDS(0, 10)
            target_cds.hsps.append(HSP(target_cds, target_domain, 100, 0, 10))
            target_cds_list.append(target_cds)

        query_cds_list = []
        for query_domain in query_domains:
            query_cds = CDS(0, 10)
            query_cds.hsps.append(HSP(query_cds, query_domain, 100, 0, 10))
            query_cds_list.append(query_cds)

        actual_results = expand_score(
            target_cds_list[target_expand_start:], query_cds_list[query_expand_start:]
        )

        self.assertEqual(expected_results, actual_results)
