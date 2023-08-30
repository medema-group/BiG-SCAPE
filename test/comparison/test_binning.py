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

        bgc_a = BGCRecord(gbk, 0, 0, 10, False, "")
        bgc_b = BGCRecord(gbk, 0, 10, 20, False, "")

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

        bgc_a = BGCRecord(gbk, 0, 0, 10, False, "")

        # b is missing GBK
        bgc_b = BGCRecord(None, 0, 0, 10, False, "")

        self.assertRaises(ValueError, BGCPair, bgc_a, bgc_b)


class TestBGCBin(TestCase):
    def test_bin_repr(self):
        """Tests whether calling str() on a bin object returns an expected string
        representation of the object
        """

        bgc_a = BGCRecord(None, 0, 0, 10, False, "")
        bgc_b = BGCRecord(None, 0, 0, 10, False, "")
        bgc_c = BGCRecord(None, 0, 0, 10, False, "")

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = BGCBin("test")

        new_bin.add_bgcs(bgc_list)

        # expected representation of the bin object
        expected_repr = "Bin 'test': 3 pairs from 3 BGCs"

        actual_repr = str(new_bin)

        self.assertEqual(expected_repr, actual_repr)

    def test_num_pairs_too_few_records(self):
        """tests if bin.num_pairs() correctly returns 0 if there is only one record in the bin"""

        gbk_a = GBK(Path("test1.gbk"), "test")
        bgc_a = BGCRecord(gbk_a, 0, 0, 10, False, "")

        new_bin = BGCBin("test")

        new_bin.add_bgcs([bgc_a])

        expected_num_pairs = 0
        actual_num_pairs = new_bin.num_pairs()

        self.assertEqual(expected_num_pairs, actual_num_pairs)

    def test_legacy_sorting(self):
        """Tests whether the legacy sorting option in bin.pairs() correctly orders the pairs"""

        gbk_a = GBK(Path("test1.gbk"), "test")
        bgc_a = BGCRecord(gbk_a, 0, 0, 10, False, "")
        gbk_b = GBK(Path("test2.gbk"), "test")
        bgc_b = BGCRecord(gbk_b, 0, 0, 10, False, "")
        gbk_c = GBK(Path("test3.gbk"), "test")
        bgc_c = BGCRecord(gbk_c, 0, 0, 10, False, "")

        # due to the order, this should generate a list of pairs as follows without legacy sort:
        # bgc_a, bgc_c
        # bgc_a, bgc_b
        # bgc_c, bgc_b
        bgc_list = [bgc_a, bgc_c, bgc_b]

        new_bin = BGCBin("test")

        new_bin.add_bgcs(bgc_list)

        # expected list should correctly sort the third entry int the list to be bgc_b, bgc_c
        expected_pair_list = [
            (bgc_a, bgc_c),
            (bgc_a, bgc_b),
            (bgc_b, bgc_c),
        ]

        actual_pair_list = [
            tuple([pair.region_a, pair.region_b])
            for pair in new_bin.pairs(None, legacy_sorting=True)
        ]

        self.assertEqual(expected_pair_list, actual_pair_list)


class TestMixComparison(TestCase):
    def test_mix_iter(self):
        """Tests whether a new mix bin can be created for comparison"""
        gbk = GBK("", "")

        bgc_a = BGCRecord(gbk, 0, 0, 10, False, "")
        bgc_a.parent_gbk = gbk

        bgc_b = BGCRecord(gbk, 0, 0, 10, False, "")
        bgc_b.parent_gbk = gbk

        bgc_c = BGCRecord(gbk, 0, 0, 10, False, "")
        bgc_c.parent_gbk = gbk

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = generate_mix(bgc_list)

        # expected representation of the bin object
        expected_pair_count = 3

        actual_pair_count = len(list(new_bin.pairs()))

        self.assertEqual(expected_pair_count, actual_pair_count)
