"""Contains tests for the calculation of comparable regions"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import GBK, Region, BGCRecord, CDS
from src.hmm import HSP
from src.comparison import BGCPair, ComparableRegion


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
            gbk_a.genes.add(cds_a)
            cds_a.hsps.add(HSP(cds_a, shared_domain, 100, 0, 30))

            cds_b = CDS(10, 90)
            gbk_b.genes.add(cds_b)
            cds_b.hsps.add(HSP(cds_b, shared_domain, 100, 0, 30))

        pair = BGCPair(gbk_a.region, gbk_b.region)

        expected_lcs = ComparableRegion(pair, 1, 4, 1, 4, False)

        actual_lcs = ComparableRegion.create_domain_lcs(pair)

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
            gbk_a.genes.add(cds_a)
            cds_a.hsps.add(HSP(cds_a, shared_domain, 100, 0, 30))

            cds_b = CDS(10, 90)
            gbk_b.genes.add(cds_b)
            cds_b.hsps.add(HSP(cds_b, shared_domain, 100, 0, 30))

        # reverse the b genes
        gbk_b.genes = gbk_b.genes[::-1]

        pair = BGCPair(gbk_a.region, gbk_b.region)

        expected_lcs = ComparableRegion(pair, 1, 4, 1, 4, True)

        actual_lcs = ComparableRegion.create_domain_lcs(pair)

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

        record = BGCRecord()
        record.nt_start = 0
        record.nt_stop = 0

        gbk = GBK(None, "test")

        non_bio_cds_1 = CDS(0, 25)
        non_bio_cds_1.gene_kind = ""
        gbk.genes.add(non_bio_cds_1)

        non_bio_cds_2 = CDS(25, 50)
        non_bio_cds_2.gene_kind = ""
        gbk.genes.add(non_bio_cds_2)

        non_bio_cds_3 = CDS(50, 75)
        non_bio_cds_3.gene_kind = ""
        gbk.genes.add(non_bio_cds_3)

        bio_cds_1 = CDS(75, 100)
        bio_cds_1.gene_kind = "biosynthetic"
        gbk.genes.add(bio_cds_1)

        record.parent_gbk = gbk

        has_biosynthetic = ComparableRegion.cds_range_contains_biosynthetic(
            record, 0, 4
        )

        self.assertTrue(has_biosynthetic)

    def test_cds_range_contains_biosyntetic_false(self):
        """tests whether the cds_range_contains_biosyntetic function returns true for a
        record in which a region contains a biosynthetic gene
        """

        record = BGCRecord()
        record.nt_start = 0
        record.nt_stop = 0

        gbk = GBK(None, "test")

        non_bio_cds_1 = CDS(0, 25)
        non_bio_cds_1.gene_kind = ""
        gbk.genes.add(non_bio_cds_1)

        non_bio_cds_2 = CDS(25, 50)
        non_bio_cds_2.gene_kind = ""
        gbk.genes.add(non_bio_cds_2)

        non_bio_cds_3 = CDS(50, 75)
        non_bio_cds_3.gene_kind = ""
        gbk.genes.add(non_bio_cds_3)

        non_bio_cds_4 = CDS(75, 100)
        non_bio_cds_4.gene_kind = ""
        gbk.genes.add(non_bio_cds_4)

        record.parent_gbk = gbk

        has_biosynthetic = ComparableRegion.cds_range_contains_biosynthetic(
            record, 0, 4
        )

        self.assertFalse(has_biosynthetic)
