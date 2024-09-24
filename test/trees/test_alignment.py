"""Contains test for GCF alignment and tree generation"""

# from python
from unittest import TestCase

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS
from big_scape.hmm import HSP, HSPAlignment
from big_scape.trees import generate_newick_tree
from big_scape.trees.newick_tree import generate_gcf_alignment, find_tree_domains


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
