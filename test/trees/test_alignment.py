"""Contains test for GCF alignment and tree generation"""

# from python
from unittest import TestCase

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS
from big_scape.hmm import HSP, HSPAlignment
from big_scape.trees import generate_newick_tree
from big_scape.trees.newick_tree import generate_gcf_alignment, find_tree_domains
from big_scape.output.legacy_output import (
    align_subrecords,
    adjust_lcs_to_family_reference,
    adjust_lcs_to_full_region,
)


class TestTrees(TestCase):
    """Contains alignment and tree generation tests"""

    def test_tree_gen_small(self):
        """Tests generated tree for families with less than three members"""

        records = [
            BGCRecord(GBK("", "", ""), 0, 0, 0, False, ""),
            BGCRecord(GBK("", "", ""), 0, 0, 0, False, ""),
        ]
        records[0]._db_id = 3
        records[1]._db_id = 64
        exemplar = 0

        expected_tree = "(3:0.0,64:0.0):0.01;"
        tree = generate_newick_tree(records, exemplar, "", "")

        self.assertEqual(tree, expected_tree)

    def test_find_tree_domains_small(self):
        """Tests whether found tree domains are correct"""
        freqs = {
            "PF1": 1,
            "PF2": 2,
            "PF3": 2,
        }

        exemplar = set(["PF1", "PF2", "PF3"])

        excepted_tree_domains = set(["PF2", "PF3"])
        actual_tree_domains = find_tree_domains(freqs, exemplar, 3)

        self.assertEqual(excepted_tree_domains, actual_tree_domains)

    def test_find_tree_domains_large(self):
        """Tests whether found tree domains are correct"""
        freqs = {
            "PF1": 5,
            "PF2": 5,
            "PF3": 7,
            "PF4": 7,
            "PF5": 6,
            "PF6": 6,
            "PF7": 4,
            "PF8": 3,
        }

        exemplar = set(["PF2", "PF3", "PF4", "PF5", "PF6", "PF7", "PF8"])

        excepted_tree_domains = set(["PF2", "PF3", "PF4", "PF5", "PF6"])
        actual_tree_domains = find_tree_domains(freqs, exemplar, 3)

        self.assertEqual(excepted_tree_domains, actual_tree_domains)

    def test_gcf_alignment(self):
        """Tests alignment of GCF HSP alignments"""
        gbk_a = GBK("", "", "")
        gbk_b = GBK("", "", "")
        cds_a = CDS(0, 20)
        cds_a.strand = 1
        cds_b = CDS(0, 20)
        cds_b.strand = 1
        hsp_a = HSP(cds_a, "PF1", 1, 0, 10)
        hsp_b = HSP(cds_b, "PF1", 1, 0, 12)
        hsp_a.alignment = HSPAlignment(hsp_a, "TEST-PF1--")
        hsp_b.alignment = HSPAlignment(hsp_b, "TEST--P-F1")
        cds_a.hsps.append(hsp_a)
        cds_b.hsps.append(hsp_b)
        gbk_a.genes.append(cds_a)
        gbk_b.genes.append(cds_b)

        records = [
            BGCRecord(gbk_a, 0, 0, 100, False, ""),
            BGCRecord(gbk_b, 1, 0, 100, False, ""),
        ]
        records[0]._db_id = 0
        records[1]._db_id = 1
        exemplar = 0
        expected_alignment = ">0\nTEST-PF1--\n>1\nTEST--P-F1\n"
        algn = generate_gcf_alignment(records, exemplar)

        self.assertEqual(algn, expected_alignment)

    def test_lcs_adjust_fwd(self):
        """Tests adjusted lcs exemplar to member not reversed"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_domain_a_start": 4,
            "lcs_domain_a_stop": 7,
            "lcs_domain_b_start": 6,
            "lcs_domain_b_stop": 9,
            "reverse": False,
        }

        expected_lcs = (4, 6, False)

        new_lcs = adjust_lcs_to_family_reference(mock_result, 0, 10, 10)

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_rev(self):
        """Tests adjusted lcs exemplar to member with reverse"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_domain_a_start": 4,
            "lcs_domain_a_stop": 7,
            "lcs_domain_b_start": 6,
            "lcs_domain_b_stop": 9,
            "reverse": True,
        }

        expected_lcs = (4, 3, True)

        new_lcs = adjust_lcs_to_family_reference(mock_result, 0, 10, 10)

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_mem2ex_fwd(self):
        """Tests adjusted lcs member to exemplar not reversed"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_domain_a_start": 6,
            "lcs_domain_a_stop": 9,
            "lcs_domain_b_start": 4,
            "lcs_domain_b_stop": 7,
            "reverse": False,
        }

        expected_lcs = (4, 6, False)

        new_lcs = adjust_lcs_to_family_reference(mock_result, 1, 10, 10)

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_mem2ex_rev(self):
        """Tests adjusted lcs member to exemplar with reverse"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_domain_a_start": 6,
            "lcs_domain_a_stop": 9,
            "lcs_domain_b_start": 4,
            "lcs_domain_b_stop": 7,
            "reverse": True,
        }

        # because the exemplar B was reversed, B is flipped back again, now
        # domain B3 corresponds to the stop in A which is corrected for exclusive start
        expected_lcs = (3, 8, True)

        new_lcs = adjust_lcs_to_family_reference(mock_result, 1, 10, 10)

        self.assertEqual(new_lcs, expected_lcs)

    def test_adjust_lcs_to_full_region_region(self):
        """Tests adjusted lcs to full regions for a region"""
        gbk_a = GBK("", "", "")
        cds_a = [CDS(100, 200), CDS(300, 550)]
        gbk_a.genes = cds_a

        region_a = BGCRecord(gbk_a, 0, 0, 600, "", "")

        gbk_b = GBK("", "", "")
        cds_b = [CDS(100, 200), CDS(300, 550), CDS(550, 700)]
        gbk_b.genes = cds_b

        region_b = BGCRecord(gbk_b, 1, 100, 700, "", "")

        a_start, b_start = (1, 1)

        expected_adjusted = (1, 1)

        actual_adjusted = adjust_lcs_to_full_region(
            a_start, b_start, region_a, region_b
        )

        self.assertEqual(expected_adjusted, actual_adjusted)

    def test_adjust_lcs_to_full_region_protocluster(self):
        """Tests adjusted lcs to full regions for a region"""
        gbk_a = GBK("", "", "")
        cds_a = [CDS(100, 200), CDS(300, 550)]
        gbk_a.genes = cds_a

        region_a = BGCRecord(gbk_a, 0, 0, 350, "", "")

        gbk_b = GBK("", "", "")
        cds_b = [CDS(100, 200), CDS(300, 550), CDS(550, 700)]
        gbk_b.genes = cds_b

        # region_b starts after the first cds
        region_b = BGCRecord(gbk_b, 1, 150, 700, "", "")

        a_start, b_start = (1, 1)

        expected_adjusted = (1, 2)

        actual_adjusted = adjust_lcs_to_full_region(
            a_start, b_start, region_a, region_b
        )

        self.assertEqual(expected_adjusted, actual_adjusted)

    def test_align_subrecords(self):
        """Tests alignment of subrecords through lcs"""
        a_domains = [
            HSP("", "PFX", 100, 0, 100),
            HSP("", "PF1", 100, 0, 100),
            HSP("", "PF2", 100, 0, 100),
            HSP("", "PF3", 100, 0, 100),
        ]

        b_domains = [
            HSP("", "PF1", 100, 0, 100),
            HSP("", "PF2", 100, 0, 100),
            HSP("", "PF3", 100, 0, 100),
            HSP("", "PFQ", 100, 0, 100),
        ]

        expected_lcs = (1, 0, False)
        actual_lcs = align_subrecords(a_domains, b_domains)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_align_subrecords_rev(self):
        """Tests alignment of subrecords through lcs"""
        a_domains = [
            HSP("", "PF3", 100, 0, 100),
            HSP("", "PF2", 100, 0, 100),
            HSP("", "PF1", 100, 0, 100),
            HSP("", "PFX", 100, 0, 100),
        ]

        b_domains = [
            HSP("", "PF1", 100, 0, 100),
            HSP("", "PF2", 100, 0, 100),
            HSP("", "PF3", 100, 0, 100),
            HSP("", "PFQ", 100, 0, 100),
        ]

        expected_lcs = (0, 2, True)
        actual_lcs = align_subrecords(a_domains, b_domains)

        self.assertEqual(expected_lcs, actual_lcs)
