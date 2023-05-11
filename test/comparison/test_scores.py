"""Contains tests for the individual and collective scores that are calculated in
BiG-SCAPE
"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import GBK, Region, CDS
from src.hmm import HSP
from src.comparison import BGCPair
from src.distances import calculate_jaccard_sets, calculate_jaccard_pair


class TestJaccard(TestCase):
    """Contains tests for the Jaccard index"""

    def test_full_overlap_set(self):
        """Tests the calculation of the Jaccard index on two overlapping sets"""
        set_a = set(["a", "b", "c", "d", "e"])
        set_b = set(["a", "b", "c", "d", "e"])

        expected_jaccard = 1.0

        actual_jaccard = calculate_jaccard_sets(set_a, set_b)

        self.assertEqual(expected_jaccard, actual_jaccard)

    def test_partial_overlap_set(self):
        """Tests the calculation of the Jaccard index on two partially overlapping sets"""
        set_a = set(["a", "b", "c", "d"])
        set_b = set(["b", "c", "d", "e"])

        expected_jaccard = 0.6

        actual_jaccard = calculate_jaccard_sets(set_a, set_b)

        self.assertAlmostEqual(expected_jaccard, actual_jaccard, 2)

    def test_no_overlap_set(self):
        """Tests the calculation of the Jaccard index on two non-overlapping sets"""
        set_a = set(["a", "b", "c"])
        set_b = set(["d", "e"])

        expected_jaccard = 0

        actual_jaccard = calculate_jaccard_sets(set_a, set_b)

        self.assertEqual(expected_jaccard, actual_jaccard)

    def test_full_overlap_pair(self):
        """Tests calculation of the jaccard index on on a BGCPair object.
        This object contains two regions of which the CDS domains overlap entirely
        """
        shared_domains = ["PF00001", "PF00002", "PF00003", "PF00004", "PF00005"]
        a_domains = [""]
        b_domains = [""]

        gbk_a = GBK("", "")
        gbk_a.region = Region(1)
        gbk_a.region.parent_gbk = gbk_a
        gbk_a.region.nt_start = 0
        gbk_a.region.nt_stop = 100

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.add(cds_a)
            cds_a.hsps.add(HSP(cds_a, a_domain, 100, 0, 30))

        gbk_b = GBK("", "")
        gbk_b.region = Region(1)
        gbk_b.region.parent_gbk = gbk_b
        gbk_b.region.nt_start = 0
        gbk_b.region.nt_stop = 100

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.add(cds_b)
            cds_b.hsps.add(HSP(cds_b, b_domain, 100, 0, 30))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            cds_b = CDS(10, 90)

            gbk_a.genes.add(cds_a)
            gbk_b.genes.add(cds_b)

            cds_a.hsps.add(HSP(cds_a, shared_domain, 100, 0, 30))
            cds_b.hsps.add(HSP(cds_b, shared_domain, 100, 0, 30))

        pair = BGCPair(gbk_a.region, gbk_b.region)

        expected_jaccard = 1

        actual_jaccard = calculate_jaccard_pair(pair)

        self.assertEqual(expected_jaccard, actual_jaccard)

    def test_partial_overlap_pair(self):
        """Tests calculation of the jaccard index on on a BGCPair object.
        This object contains two regions of which the CDS domains overlap partially
        """
        shared_domains = [
            "PF00002",
            "PF00003",
            "PF00004",
        ]
        a_domains = [
            "PF00001",
        ]
        b_domains = [
            "PF00005",
        ]

        gbk_a = GBK("", "")
        gbk_a.region = Region(1)
        gbk_a.region.parent_gbk = gbk_a
        gbk_a.region.nt_start = 0
        gbk_a.region.nt_stop = 100

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.add(cds_a)
            cds_a.hsps.add(HSP(cds_a, a_domain, 100, 0, 30))

        gbk_b = GBK("", "")
        gbk_b.region = Region(1)
        gbk_b.region.parent_gbk = gbk_b
        gbk_b.region.nt_start = 0
        gbk_b.region.nt_stop = 100

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.add(cds_b)
            cds_b.hsps.add(HSP(cds_b, b_domain, 100, 0, 30))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            cds_b = CDS(10, 90)

            gbk_a.genes.add(cds_a)
            gbk_b.genes.add(cds_b)

            cds_a.hsps.add(HSP(cds_a, shared_domain, 100, 0, 30))
            cds_b.hsps.add(HSP(cds_b, shared_domain, 100, 0, 30))

        pair = BGCPair(gbk_a.region, gbk_b.region)

        expected_jaccard = 0.6

        actual_jaccard = calculate_jaccard_pair(pair)

        self.assertAlmostEqual(expected_jaccard, actual_jaccard, 2)

    def test_no_overlap_pair(self):
        """Tests calculation of the jaccard index on on a BGCPair object.
        This object contains two regions of which the CDS domains do not overlap
        """
        shared_domains = []
        a_domains = [
            "PF00001",
            "PF00002",
            "PF00003",
        ]
        b_domains = [
            "PF00005",
            "PF00004",
        ]

        gbk_a = GBK("", "")
        gbk_a.region = Region(1)
        gbk_a.region.parent_gbk = gbk_a
        gbk_a.region.nt_start = 0
        gbk_a.region.nt_stop = 100

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.add(cds_a)
            cds_a.hsps.add(HSP(cds_a, a_domain, 100, 0, 30))

        gbk_b = GBK("", "")
        gbk_b.region = Region(1)
        gbk_b.region.parent_gbk = gbk_b
        gbk_b.region.nt_start = 0
        gbk_b.region.nt_stop = 100

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.add(cds_b)
            cds_b.hsps.add(HSP(cds_b, b_domain, 100, 0, 30))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            cds_b = CDS(10, 90)

            gbk_a.genes.add(cds_a)
            gbk_b.genes.add(cds_b)

            cds_a.hsps.add(HSP(cds_a, shared_domain, 100, 0, 30))
            cds_b.hsps.add(HSP(cds_b, shared_domain, 100, 0, 30))

        pair = BGCPair(gbk_a.region, gbk_b.region)

        expected_jaccard = 0

        actual_jaccard = calculate_jaccard_pair(pair)

        self.assertEqual(expected_jaccard, actual_jaccard, 2)


class TestAdjacency(TestCase):
    """Contains tests for the Jaccard index"""


class TestDSS(TestCase):
    """Contains tests for the Jaccard index"""
