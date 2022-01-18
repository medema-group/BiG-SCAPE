from unittest import TestCase

from src import big_scape

class TestRunBase(TestCase):
    """Class containing tests for the run object, which tracks run details and keeps timing
    """
    test_run: big_scape.Run

    def setUp(self):
    
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

        # results in float, so asset almost equal
        self.assertAlmostEqual(jaccard_index, test_jaccard_index)

    def test_get_distance_default(self):
        """Tests the distance metrics returned by calc_distance_lcs on a pair of mock clusters
        
        This is the default version of this test, and is equivalent to a run of BiG-SCAPE with
        no extra command line arguments that have an effect in calc_distance_lcs
        """
        # test run details
        # we have to fill this in ourselves until there is a proper mocking mechanic in place for
        # this object
        test_run = big_scape.Run()


        distance = big_scape.scores.calc_distance_lcs(run, cluster_info_a, cluster_info_b,
                                                      weights, bgcs, domain_count_gene
                                                      bgc_gene_orientation, bgc_info,
                                                      aligned_domain_sequences)
