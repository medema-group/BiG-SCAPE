"""Contains test for GCF alignment and tree generation"""

# from python
from unittest import TestCase

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS
from big_scape.hmm import HSP, HSPAlignment
from big_scape.trees import generate_newick_tree
from big_scape.trees.newick_tree import generate_gcf_alignment
from big_scape.output.legacy_output import adjust_lcs_to_all_genes


class TestTrees(TestCase):
    """Contains alignment and tree generation tests"""

    def test_tree_gen_small(self):
        """Tests generated tree for families with less than three members"""

        records = [
            BGCRecord(GBK("", "", ""), 0, 0, 0, False, ""),
            BGCRecord(GBK("", "", ""), 0, 0, 0, False, ""),
        ]
        exemplar = 0
        mock_family = [0, 1]

        expected_tree = "(0:0.0,1:0.0):0.01;"
        tree = generate_newick_tree(records, exemplar, mock_family)

        self.assertEqual(tree, expected_tree)

    def test_gcf_alignment(self):
        """Tests alignment of GCF HSP alignments"""
        gbk_a = GBK("", "", "")
        gbk_b = GBK("", "", "")
        cds_a = CDS(0, 20)
        cds_b = CDS(0, 20)
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
        exemplar = 0
        family_members = [0, 1]
        expected_alignment = ">0\nTEST-PF1--\n>1\nTEST--P-F1\n"
        algn = generate_gcf_alignment(records, exemplar, family_members)

        self.assertEqual(algn, expected_alignment)

    def test_lcs_adjust_ex2mem(self):
        """Tests adjusted lcs exemplar to member not reversed"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 4,
            "lcs_a_stop": 7,
            "lcs_b_start": 6,
            "lcs_b_stop": 9,
            "reverse": False,
        }
        dom_to_gene = {0: {4: 5, 7: 9}, 1: {6: 8}}
        dom_count = {}

        expected_lcs = (5, 8, False)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 0, 1, GBK("", "", ""), GBK("", "", ""), dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_ex2mem_reverse(self):
        """Tests adjusted lcs exemplar to member with reverse"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 2,
            "lcs_a_stop": 3,
            "lcs_b_start": 1,
            "lcs_b_stop": 2,
            "reverse": True,
        }
        dom_to_gene = {0: {0: 0, 1: 2, 2: 4}, 1: {0: 1, 1: 2, 2: 3}}
        dom_count = {0: [1, 1, 1], 1: [1, 1, 1]}

        expected_lcs = (4, 2, True)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 0, 1, GBK("", "", ""), GBK("", "", ""), dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_mem2ex(self):
        """Tests adjusted lcs member to exemplar not reversed"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 2,
            "lcs_a_stop": 3,
            "lcs_b_start": 1,
            "lcs_b_stop": 2,
            "reverse": False,
        }
        dom_to_gene = {0: {0: 0, 1: 2, 2: 4}, 1: {0: 1, 1: 2, 2: 3}}
        dom_count = {0: [1, 1, 1], 1: [1, 1, 1]}

        expected_lcs = (2, 4, False)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 1, 0, GBK("", "", ""), GBK("", "", ""), dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_mem2ex_reverse(self):
        """Tests adjusted lcs member to exemplar with reverse"""
        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 2,
            "lcs_a_stop": 3,
            "lcs_b_start": 1,
            "lcs_b_stop": 2,
            "reverse": True,
        }
        dom_to_gene = {0: {0: 0, 1: 2, 2: 4}, 1: {0: 1, 1: 2, 2: 3}}
        dom_count = {0: [1, 1, 1], 1: [1, 1, 1]}

        expected_lcs = (2, 4, True)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 1, 0, GBK("", "", ""), GBK("", "", ""), dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_zero_length_not_reversed(self):
        """Tests adjusted lcs for lcs with length 0"""

        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 0,
            "lcs_a_stop": 0,
            "lcs_b_start": 0,
            "lcs_b_stop": 0,
            "reverse": False,
        }
        dom_to_gene = {0: {0: 0, 1: 2, 2: 4}, 1: {0: 1, 1: 2, 2: 3}}
        dom_count = {0: [1, 2, 1], 1: [1, 1, 2]}

        max_dom_cds = CDS(20, 30)
        max_dom_cds.strand = 1
        fill_cds = CDS(10, 20)
        fill_cds.strand = -1

        exempl_gbk = GBK("", "", "")
        exempl_gbk.genes = [fill_cds, fill_cds, max_dom_cds, fill_cds, fill_cds]

        mem_gbk = GBK("", "", "")
        mem_gbk.genes = [CDS(0, 1), CDS(2, 3), CDS(4, 5), max_dom_cds]

        expected_lcs = (2, 3, False)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 0, 1, exempl_gbk, mem_gbk, dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)

    def test_lcs_adjust_zero_length_reversed(self):
        """Tests adjusted lcs for lcs with length 0"""

        mock_result = {
            "record_a_id": 0,
            "record_b_id": 1,
            "lcs_a_start": 0,
            "lcs_a_stop": 0,
            "lcs_b_start": 0,
            "lcs_b_stop": 0,
            "reverse": False,
        }
        dom_to_gene = {0: {0: 0, 1: 2, 2: 4}, 1: {0: 1, 1: 2, 2: 3}}
        dom_count = {0: [1, 2, 1], 1: [1, 1, 2]}

        max_dom_cds = CDS(0, 0)
        max_dom_cds.strand = 1
        fill_cds = CDS(0, 0)
        fill_cds.strand = -1

        exempl_gbk = GBK("", "", "")
        exempl_gbk.genes = [fill_cds, fill_cds, max_dom_cds, fill_cds, fill_cds]

        max_dom_cds_rev = CDS(0, 0)
        max_dom_cds_rev.strand = -1
        mem_gbk = GBK("", "", "")
        mem_gbk.genes = [CDS(0, 0), CDS(0, 0), CDS(0, 0), max_dom_cds_rev]

        expected_lcs = (2, 1, True)

        new_lcs = adjust_lcs_to_all_genes(
            mock_result, 0, 1, exempl_gbk, mem_gbk, dom_to_gene, dom_count
        )

        self.assertEqual(new_lcs, expected_lcs)
