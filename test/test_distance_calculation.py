from unittest import TestCase

from src import big_scape

class TestRunBase(TestCase):
    """Class containing tests for the run object, which tracks run details and keeps timing
    """
    test_run: big_scape.Run

    def setUp(self):
        # test run object for various parameters
        self.test_run = big_scape.Run()
    
    def test_jaccard_index(self):
        """Tests the jaccard index caluculation """
        # test sets. there are three elements which are common, and 2 elements which are not in
        # common
        test_set_a = set([1, 2, 3, 4, 5])
        test_set_b = set([0, 2, 3, 4, 6])

        # the jaccard index is defined as the number of intersecting elements divided by the number 
        # of elements in the union of two sets (A ∩ B / B ∪ B)
        intersect = test_set_a & test_set_b
        union = test_set_a | test_set_b

        intersect_len = len(intersect)
        union_len = len(union)

        test_jaccard_index = intersect_len / union_len

        jaccard_index = big_scape.scores.calc_jaccard(intersect, union)

        self.assertAlmostEqual(jaccard_index, test_jaccard_index)
