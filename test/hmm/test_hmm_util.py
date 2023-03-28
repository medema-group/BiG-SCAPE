"""Contains tests for the utility functions in the hmm module"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import CDS
from src.hmm import HSP
from src.hmm.hmm import has_overlap, len_overlap


class TestHMMUtil(TestCase):
    """Contains tests to cover any utility functions in the hmm module"""

    def test_no_overlap_left(self):
        """Tests the no_overlap function where a is left of b"""
        cds_a = CDS(0, 50)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = False

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_no_overlap_right(self):
        """Tests the no_overlap function where a is right of b"""
        cds_a = CDS(100, 150)
        cds_b = CDS(0, 50)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = False

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_full(self):
        """Tests the no_overlap function for two overlapping regions

        both regions overlap perfectly"""
        cds_a = CDS(100, 150)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_partial_left(self):
        """Tests the no_overlap function for two overlapping regions

        a overlaps b on b left side"""
        cds_a = CDS(0, 120)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_partial_right(self):
        """Tests the no_overlap function for two overlapping regions

        a overlaps b on b right side"""
        cds_a = CDS(130, 250)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_full(self):
        """Tests the len_overlap function for two fully overlapping regions"""
        cds_a = CDS(0, 100)
        cds_b = CDS(0, 100)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 100

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_left(self):
        """Tests the len_overlap function for two overlapping regions with 50 bp overlap

        a overlaps b on b left side"""
        cds_a = CDS(0, 100)
        cds_b = CDS(50, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 50

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_right(self):
        """Tests the len_overlap function for two overlapping regions with 50 bp overlap

        a overlaps b on b right side"""
        cds_a = CDS(100, 200)
        cds_b = CDS(50, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 50

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_no_overlap_left(self):
        """Tests the len_overlap function for two non-overlapping regions

        a < b"""
        cds_a = CDS(0, 100)
        cds_b = CDS(100, 200)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 0

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_no_overlap_right(self):
        """Tests the len_overlap function for two non-overlapping regions

        a > b"""
        cds_a = CDS(200, 300)
        cds_b = CDS(100, 200)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 0

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)
