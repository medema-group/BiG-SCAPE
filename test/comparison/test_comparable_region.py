"""Contains tests for the calculation of comparable regions"""

# from python
from unittest import TestCase

# from other modules
from big_scape.genbank import GBK, Region, CDS
from big_scape.hmm import HSP
from big_scape.comparison import RecordPair, ComparableRegion


class TestComparibleRegions(TestCase):
    """Contains tests for calulation of comparable regions between two BGCs"""

    def test_repr_fwd(self):
        """Tests the comparable region object string representation when B is not
        reversed
        """
        comparable_region = ComparableRegion(0, 5, 3, 8, 0, 0, 0, 0, False)

        expected_repr = "Comparable region: A 0-5, B 3-8, B is not reversed"

        actual_repr = str(comparable_region)

        self.assertEqual(expected_repr, actual_repr)

    def test_repr_rev(self):
        """Tests the comparable region object string representation when B is not
        reversed
        """
        comparable_region = ComparableRegion(0, 5, 3, 8, 0, 0, 0, 0, True)

        expected_repr = "Comparable region: A 0-5, B 3-8, B is reversed"

        actual_repr = str(comparable_region)

        self.assertEqual(expected_repr, actual_repr)

    def test_cds_range_contains_biosyntetic_true(self):
        """tests whether the cds_range_contains_biosyntetic function returns true for a
        record in which a region contains a biosynthetic gene
        """

        gbk = GBK(None, "test", "test")

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

        gbk = GBK(None, "test", "test")

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

        gbk_a = GBK("", "", "")
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")

        for a_domain in a_domains:
            cds_a = CDS(10, 90)
            gbk_a.genes.append(cds_a)
            env_start = len(gbk_a.genes) * 10
            env_stop = env_start + 10
            cds_a.hsps.append(HSP(cds_a, a_domain, 100, env_start, env_stop))

        gbk_b = GBK("", "", "")
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

        actual_dicts = pair.get_domain_dicts()

        self.assertEqual(expected_dicts, actual_dicts)

    def test_inflate(self):
        """Tests inflate function of comparable region"""
        # Two gbks with five cds
        # only middle gbk A cds contain domains
        gbk_a = GBK("", "", "test")
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")
        cds_a = [CDS(10, 20), CDS(30, 40), CDS(50, 60), CDS(70, 80), CDS(80, 90)]
        # assign orf numbers
        for idx, cds in enumerate(cds_a):
            cds.orf_num = idx + 1
        cds_a[1].hsps.append(HSP(cds_a[1], "PF1", 100, 0, 10))
        cds_a[2].hsps.append(HSP(cds_a[2], "PF1", 100, 0, 10))
        cds_a[3].hsps.append(HSP(cds_a[3], "PF1", 100, 0, 10))
        gbk_a.genes = cds_a

        # repeat for gbk B, only start and end cds of gbk_b contain domain
        gbk_b = GBK("", "", "test")
        gbk_b.region = Region(gbk_b, 0, 0, 100, False, "")
        cds_b = [CDS(10, 20), CDS(30, 40), CDS(50, 60), CDS(70, 80), CDS(80, 90)]
        for idx, cds in enumerate(cds_b):
            cds.orf_num = idx + 1
        cds_b[0].hsps.append(HSP(cds_b[0], "PF1", 100, 0, 10))
        cds_b[4].hsps.append(HSP(cds_b[4], "PF1", 100, 0, 10))
        gbk_b.genes = cds_b

        pair = RecordPair(gbk_a.region, gbk_b.region)

        # set mock lcs and extend bounds (stops are exclusive)
        pair.comparable_region.lcs_a_start = 0
        pair.comparable_region.lcs_a_stop = 2
        pair.comparable_region.a_start = 1
        pair.comparable_region.a_stop = 3

        pair.comparable_region.lcs_b_start = 0
        pair.comparable_region.lcs_b_stop = 2
        pair.comparable_region.b_start = 0
        pair.comparable_region.b_stop = 2
        pair.comparable_region.reverse = False

        expected_inflate = (1, 3, 2, 4, 0, 5, 0, 5)

        pair.comparable_region.inflate(pair)

        actual_inflate = (
            pair.comparable_region.lcs_a_start,
            pair.comparable_region.lcs_a_stop,
            pair.comparable_region.a_start,
            pair.comparable_region.a_stop,
            pair.comparable_region.lcs_b_start,
            pair.comparable_region.lcs_b_stop,
            pair.comparable_region.b_start,
            pair.comparable_region.b_stop,
        )

        self.assertEqual(expected_inflate, actual_inflate)

    def test_inflate_b_reverse(self):
        """Tests infate_b when reverse is true"""
        gbk_a = GBK("", "", "test")
        gbk_a.region = Region(gbk_a, 0, 0, 100, False, "")

        gbk_b = GBK("", "", "test")
        gbk_b.region = Region(gbk_b, 0, 0, 100, False, "")
        cds_b = [CDS(10, 20), CDS(30, 40), CDS(50, 60), CDS(70, 80), CDS(80, 90)]
        for idx, cds in enumerate(cds_b):
            cds.orf_num = idx + 1
        cds_b[1].hsps.append(HSP(cds_b[1], "PF1", 100, 0, 10))
        cds_b[2].hsps.append(HSP(cds_b[2], "PF1", 100, 0, 10))
        gbk_b.genes = cds_b

        pair = RecordPair(gbk_a.region, gbk_b.region)

        # only interested in B, reverse=True (stops are exclusive)
        pair.comparable_region.lcs_b_start = 0
        pair.comparable_region.lcs_b_stop = 2
        pair.comparable_region.b_start = 0
        pair.comparable_region.b_stop = 2
        pair.comparable_region.reverse = True

        expected_inflate = (2, 4, 2, 4)

        pair.comparable_region.inflate_b(pair.record_b)

        actual_inflate = (
            pair.comparable_region.lcs_b_start,
            pair.comparable_region.lcs_b_stop,
            pair.comparable_region.b_start,
            pair.comparable_region.b_stop,
        )

        self.assertEqual(expected_inflate, actual_inflate)
