"""Contains tests for the pair class and its methods"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from src.genbank import GBK, BGCRecord
from src.comparison import BGCBin, BGCPair
from src.comparison import generate_mix


class TestBGCPair(TestCase):
    """Contains tests to test pair objects used in binning"""

    def test_pair_repr(self):
        """Tests whether calling str() on a bin object returns an expected string
        representation of the object
        """
        gbk = GBK(Path("test"), "test")

        bgc_a = BGCRecord()
        bgc_a.parent_gbk = gbk
        bgc_a.nt_start = 0
        bgc_a.nt_stop = 10

        bgc_b = BGCRecord()
        bgc_b.parent_gbk = gbk
        bgc_b.nt_start = 10
        bgc_b.nt_stop = 20

        pair = BGCPair(bgc_a, bgc_b)

        expected_repr = (
            "Pair GBK test, 0 genes Record (superclass) 0-10 - GBK test, 0 genes "
            "Record (superclass) 10-20"
        )

        actual_repr = str(pair)

        self.assertEqual(expected_repr, actual_repr)

    def test_pair_no_parent_gbk(self):
        """Tests whether initialization of a BGC pair where one of the BGCs does not
        have a parent GBK correctly throws a ValueError
        """
        gbk = GBK("", "")

        bgc_a = BGCRecord()
        bgc_a.parent_gbk = gbk

        # b is missing GBK
        bgc_b = BGCRecord()

        self.assertRaises(ValueError, BGCPair, bgc_a, bgc_b)


class TestBGCBin(TestCase):
    def test_bin_repr(self):
        """Tests whether calling str() on a bin object returns an expected string
        representation of the object
        """

        bgc_a = BGCRecord()
        bgc_b = BGCRecord()
        bgc_c = BGCRecord()

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = BGCBin("test")

        new_bin.add_bgcs(bgc_list)

        # expected representation of the bin object
        expected_repr = "Bin 'test': 3 pairs from 3 BGCs"

        actual_repr = str(new_bin)

        self.assertEqual(expected_repr, actual_repr)


class TestMixComparison(TestCase):
    def test_mix_iter(self):
        """Tests whether a new mix bin can be created for comparison"""
        gbk = GBK("", "")

        bgc_a = BGCRecord()
        bgc_a.parent_gbk = gbk

        bgc_b = BGCRecord()
        bgc_b.parent_gbk = gbk

        bgc_c = BGCRecord()
        bgc_c.parent_gbk = gbk

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = generate_mix(bgc_list)

        # expected representation of the bin object
        expected_pair_count = 3

        actual_pair_count = len(list(new_bin.pairs()))

        self.assertEqual(expected_pair_count, actual_pair_count)
