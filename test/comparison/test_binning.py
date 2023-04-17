"""Contains tests for the pair class and its methods"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import BGCRecord
from src.comparison import generate_mix


class TestComparisonBin(TestCase):
    def test_create_mix(self):
        """Tests whether a new mix bin can be created for comparison"""

        bgc_a = BGCRecord()
        bgc_b = BGCRecord()
        bgc_c = BGCRecord()

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = generate_mix(bgc_list)

        # expected representation of the bin object
        expected_repr = "Bin 'mix': 3 pairs from 3 BGCs"

        actual_repr = str(new_bin)

        self.assertEqual(expected_repr, actual_repr)
