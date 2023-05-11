"""Contains tests for the individual and collective scores that are calculated in
BiG-SCAPE
"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import GBK, Region, CDS
from src.hmm import HSP
from src.comparison import BGCPair
from src.distances import (
    calc_jaccard_sets,
    calc_jaccard_pair,
    calc_ai_lists,
    calc_ai_pair,
)


class TestJaccard(TestCase):
    """Contains tests for the Jaccard index"""

    def test_full_overlap_set(self):
        """Tests the calculation of the Jaccard index on two overlapping sets"""
        set_a = set(["a", "b", "c", "d", "e"])
        set_b = set(["a", "b", "c", "d", "e"])

        expected_jaccard = 1.0

        actual_jaccard = calc_jaccard_sets(set_a, set_b)

        self.assertEqual(expected_jaccard, actual_jaccard)

    def test_partial_overlap_set(self):
        """Tests the calculation of the Jaccard index on two partially overlapping sets"""
        set_a = set(["a", "b", "c", "d"])
        set_b = set(["b", "c", "d", "e"])

        expected_jaccard = 0.6

        actual_jaccard = calc_jaccard_sets(set_a, set_b)

        self.assertAlmostEqual(expected_jaccard, actual_jaccard, 2)

    def test_no_overlap_set(self):
        """Tests the calculation of the Jaccard index on two non-overlapping sets"""
        set_a = set(["a", "b", "c"])
        set_b = set(["d", "e"])

        expected_jaccard = 0

        actual_jaccard = calc_jaccard_sets(set_a, set_b)

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

        actual_jaccard = calc_jaccard_pair(pair)

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

        actual_jaccard = calc_jaccard_pair(pair)

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

        actual_jaccard = calc_jaccard_pair(pair)

        self.assertEqual(expected_jaccard, actual_jaccard, 2)


class TestAdjacency(TestCase):
    """Contains tests for the Adjacency index"""

    def test_full_adjacent_lists(self):
        """Tests the calculation of the Adjacency index on two lists with the same
        adjacent items
        """
        list_a = ["a", "b", "c", "d", "e"]
        list_b = ["a", "b", "c", "d", "e"]

        # a sets: ab, bc, cd, de
        # b sets: ab, bc, cd, de
        # intersect = ab, bc, cd, de (4)
        # union = ab, bc, cd, de (4)
        # adjacency = 4 / 4 = 1.0

        expected_adjacency = 1.0

        actual_adjacency = calc_ai_lists(list_a, list_b)

        self.assertAlmostEqual(expected_adjacency, actual_adjacency, 2)

    def test_partial_adjacent_lists(self):
        """Tests the calculation of the adjacency index on two lists with items that are
        partialy adjacent"""
        # ab adjacent for both, cd adjacent for both, but bc not
        list_a = ["a", "b", "c", "d"]
        list_b = ["c", "d", "a", "b"]

        # a sets: ab, bc, cd
        # b sets: cd, da, ab
        # intersect = ab & cd (2)
        # union = ab, bc, da, cd (4)
        # adjacency = 2 / 4 = 0.5

        expected_adjacency = 0.5

        actual_adjacency = calc_ai_lists(list_a, list_b)

        self.assertAlmostEqual(expected_adjacency, actual_adjacency, 2)

    def test_no_adjacent_list(self):
        """Tests the calculation of the Jaccard index on two non-overlapping sets"""
        list_a = ["a", "b", "c"]
        list_b = ["d", "e", "f"]

        # a sets: ab, bc
        # b sets: de, ef
        # intersect = None (0)
        # union = ab, bc, de, ef (4)
        # adjacency = 0 / 4 = 0

        expected_adjacency = 0.0

        actual_adjacency = calc_ai_lists(list_a, list_b)

        self.assertAlmostEqual(expected_adjacency, actual_adjacency, 2)

    def test_all_adjacent_pair(self):
        """Tests calculation of the jaccard index on on a BGCPair object.
        This object contains two regions of which the CDS domains overlap partially
        """
        shared_domains = [
            "PF00001",
            "PF00002",
            "PF00003",
        ]
        a_domains = []
        b_domains = []

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

        expected_ai = 1.0

        actual_ai = calc_ai_pair(pair)

        self.assertAlmostEqual(expected_ai, actual_ai, 2)


class TestDSS(TestCase):
    """Contains tests for the Domain similarity index"""
