"""Contains tests for the calculation of comparable regions"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import GBK, Region, CDS
from src.hmm import HSP
from src.comparison import RecordPair, ComparableRegion


class TestComparibleRegions(TestCase):
    """Contains tests for calulation of comparable regions between two BGCs"""

    def test_get_dom_list_lcs(self):
        """Tests whether get_dom_list_lcs can find the longest common substring of
        domains between a pair of BGCs"""
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
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            cds_a.hsps.append(HSP(cds_a, a_domain, 100, 0, 30))

        gbk_b = GBK("", "")
        gbk_b.region = Region(gbk_b, 0, 0, 100, False, "")

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            cds_b.hsps.append(HSP(cds_b, b_domain, 100, 0, 30))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            cds_a.hsps.append(HSP(cds_a, shared_domain, 100, 0, 30))

            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            cds_b.hsps.append(HSP(cds_b, shared_domain, 100, 0, 30))

        pair = RecordPair(gbk_a.region, gbk_b.region)

        expected_lcs = ComparableRegion(pair, 1, 4, 1, 4, False)

        pair.comparable_region.find_lcs()

        actual_lcs = pair.comparable_region

        self.assertEqual(expected_lcs, actual_lcs)

    def test_get_dom_list_lcs_reverse(self):
        """Tests whether get_dom_list_lcs can find the longest common substring of
        domains between a pair of BGCs where one side has a reversed domain string"""
        a_domains = [
            "PF00001",
        ]
        b_domains = [
            "PF00005",
        ]
        shared_domains = [
            "PF00002",
            "PF00003",
            "PF00004",
        ]

        gbk_a = GBK("", "")
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            cds_a.hsps.append(HSP(cds_a, a_domain, 100, 0, 30))

        gbk_b = GBK("", "")
        gbk_b.region = Region(gbk_b, 0, 0, 100, False, "")

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            cds_b.hsps.append(HSP(cds_b, b_domain, 100, 0, 30))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            cds_a.hsps.append(HSP(cds_a, shared_domain, 100, 0, 30))

            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            cds_b.hsps.append(HSP(cds_b, shared_domain, 100, 0, 30))

        # reverse the b genes
        gbk_b.genes = gbk_b.genes[::-1]

        pair = RecordPair(gbk_a.region, gbk_b.region)

        expected_lcs = ComparableRegion(pair, 1, 4, 1, 4, True)

        pair.comparable_region.find_lcs()

        actual_lcs = pair.comparable_region

        self.assertEqual(expected_lcs, actual_lcs)

    def test_repr_fwd(self):
        """Tests the comparable region object string representation when B is not
        reversed
        """
        comparable_region = ComparableRegion(None, 0, 5, 3, 8, False)

        expected_repr = "Comparable region: A 0-5, B 3-8, B is not reversed"

        actual_repr = str(comparable_region)

        self.assertEqual(expected_repr, actual_repr)

    def test_repr_rev(self):
        """Tests the comparable region object string representation when B is not
        reversed
        """
        comparable_region = ComparableRegion(None, 0, 5, 3, 8, True)

        expected_repr = "Comparable region: A 0-5, B 3-8, B is reversed"

        actual_repr = str(comparable_region)

        self.assertEqual(expected_repr, actual_repr)

    def test_cds_range_contains_biosyntetic_true(self):
        """tests whether the cds_range_contains_biosyntetic function returns true for a
        record in which a region contains a biosynthetic gene
        """

        gbk = GBK(None, "test")

        record = Region(gbk, 0, 0, 0, False, "")

        non_bio_cds_1 = CDS(0, 25)
        non_bio_cds_1.hsps.append(HSP(non_bio_cds_1, "test", 100.0, 0, 25))
        non_bio_cds_1.gene_kind = ""
        gbk.genes.append(non_bio_cds_1)

        non_bio_cds_2 = CDS(25, 50)
        non_bio_cds_2.hsps.append(HSP(non_bio_cds_2, "test", 100.0, 0, 25))
        non_bio_cds_2.gene_kind = ""
        gbk.genes.append(non_bio_cds_2)

        non_bio_cds_3 = CDS(50, 75)
        non_bio_cds_3.hsps.append(HSP(non_bio_cds_3, "test", 100.0, 0, 25))
        non_bio_cds_3.gene_kind = ""
        gbk.genes.append(non_bio_cds_3)

        bio_cds_1 = CDS(75, 100)
        bio_cds_1.hsps.append(HSP(bio_cds_1, "test", 100.0, 0, 25))
        bio_cds_1.gene_kind = "biosynthetic"
        gbk.genes.append(bio_cds_1)

        record.parent_gbk = gbk

        has_biosynthetic = ComparableRegion.cds_range_contains_biosynthetic(
            record, 0, 4
        )

        self.assertTrue(has_biosynthetic)

    def test_cds_range_contains_biosyntetic_false(self):
        """tests whether the cds_range_contains_biosyntetic function returns true for a
        record in which a region contains a biosynthetic gene
        """

        gbk = GBK(None, "test")

        record = Region(gbk, 0, 0, 0, False, "")

        non_bio_cds_1 = CDS(0, 25)
        non_bio_cds_1.gene_kind = ""
        gbk.genes.append(non_bio_cds_1)

        non_bio_cds_2 = CDS(25, 50)
        non_bio_cds_2.gene_kind = ""
        gbk.genes.append(non_bio_cds_2)

        non_bio_cds_3 = CDS(50, 75)
        non_bio_cds_3.gene_kind = ""
        gbk.genes.append(non_bio_cds_3)

        non_bio_cds_4 = CDS(75, 100)
        non_bio_cds_4.gene_kind = ""
        gbk.genes.append(non_bio_cds_4)

        record.parent_gbk = gbk

        has_biosynthetic = ComparableRegion.cds_range_contains_biosynthetic(
            record, 0, 4
        )

        self.assertFalse(has_biosynthetic)

    def test_get_domain_dicts(self):
        """Tests whether get_dom_dict will return the correct dictionaries of domains"""
        shared_domains = [
            "PF00002",
            "PF00003",
            "PF00004",
        ]
        a_domains = [
            "PF00001",
            "PF00002",
        ]
        b_domains = [
            "PF00005",
            "PF00004",
        ]

        gbk_a = GBK("", "")
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            env_start = len(gbk_a.genes) * 10
            env_stop = env_start + 10
            cds_a.hsps.append(HSP(cds_a, a_domain, 100, env_start, env_stop))

        gbk_b = GBK("", "")
        gbk_b.region = Region(gbk_b, 0, 0, 100, False, "")

        for b_domain in b_domains:
            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            env_start = len(gbk_b.genes) * 10
            env_stop = env_start + 10
            cds_b.hsps.append(HSP(cds_b, b_domain, 100, env_start, env_stop))

        for shared_domain in shared_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            env_start = len(gbk_a.genes) * 10
            env_stop = env_start + 10
            cds_a.hsps.append(HSP(cds_a, shared_domain, 100, env_start, env_stop))

            cds_b = CDS(10, 90)
            gbk_b.genes.append(cds_b)
            env_start = len(gbk_b.genes) * 10
            env_stop = env_start + 10
            cds_b.hsps.append(HSP(cds_b, shared_domain, 100, env_start, env_stop))

        pair = RecordPair(gbk_a.region, gbk_b.region)

        expected_dicts = (
            {
                "PF00001": [0],
                "PF00002": [1, 2],
                "PF00003": [3],
                "PF00004": [4],
            },
            {
                "PF00004": [1, 4],
                "PF00003": [3],
                "PF00002": [2],
                "PF00005": [0],
            },
        )

        actual_dicts = pair.comparable_region.get_domain_dicts()

        self.assertEqual(expected_dicts, actual_dicts)
