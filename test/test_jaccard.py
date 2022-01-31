from unittest import TestCase
from src.big_scape.bgc_info import BgcInfo

from src import big_scape

class TestJaccard(TestCase):
    """Test class for jaccard distance.

    This class constructs two couples of bgcs
    one couple has overlapping domains, the other does not

    Jaccard distance should be dissimilar between these two couples
    """
    test_run: big_scape.Run

    similar_bgc_a: BgcInfo
    similar_bgc_b: BgcInfo

    distinct_bgc_a: BgcInfo
    distinct_bgc_b: BgcInfo

    def setUp(self):
        self.similar_bgc_a = BgcInfo("test_similar_bgc_a")
        self.similar_bgc_b = BgcInfo("test_similar_bgc_b")

        self.similar_bgc_a.ordered_domain_list = []
        self.similar_bgc_b.ordered_domain_list = []

        # we'll add 10 domains to these genes
        for i in range(10):
            dom_name = "DOM_" + str(i)
            # bottom 8 are in a
            if i < 8:
                self.similar_bgc_a.ordered_domain_list.append(dom_name)
            # top 8 are in b
            if i > 1:
                self.similar_bgc_b.ordered_domain_list.append(dom_name)
        
        # for dissimilar we'll do the same but with less overlap
        self.distinct_bgc_a = BgcInfo("test_distinct_bgc_a")
        self.distinct_bgc_b = BgcInfo("test_distinct_bgc_b")

        self.distinct_bgc_a.ordered_domain_list = []
        self.distinct_bgc_b.ordered_domain_list = []

        # we'll add 10 domains to these genes
        for i in range(15):
            dom_name = "DOM_" + str(i)
            # even are in a
            if i % 2 == 0:
                self.distinct_bgc_a.ordered_domain_list.append(dom_name)
            # odd are in b
            else:
                self.distinct_bgc_b.ordered_domain_list.append(dom_name)

    
    def test_jaccard_index_similar(self):
        """Tests the jaccard index caluculation """
        # test sets. there are three elements which are common, and 2 elements which are not in
        # common
        test_set_a = set(self.similar_bgc_a.ordered_domain_list)
        test_set_b = set(self.similar_bgc_b.ordered_domain_list)

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


    def test_jaccard_index_distinct(self):
        """Tests the jaccard index caluculation """
        # test sets. there are three elements which are common, and 2 elements which are not in
        # common
        test_set_a = set(self.distinct_bgc_a.ordered_domain_list)
        test_set_b = set(self.distinct_bgc_b.ordered_domain_list)

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
    

    def test_jaccard_index_compare(self):
        test_set_a_similar = set(self.similar_bgc_a.ordered_domain_list)
        test_set_b_similar = set(self.similar_bgc_b.ordered_domain_list)
        
        test_set_a_distinct = set(self.distinct_bgc_a.ordered_domain_list)
        test_set_b_distinct = set(self.distinct_bgc_b.ordered_domain_list)

        similar_intersect = test_set_a_similar & test_set_b_similar
        similar_union = test_set_a_similar | test_set_b_similar
        similar_jaccard = big_scape.scores.calc_jaccard(similar_intersect, similar_union)

        distinct_intersect = test_set_a_distinct & test_set_b_distinct
        distinct_union = test_set_a_distinct | test_set_b_distinct
        distinct_jaccard = big_scape.scores.calc_jaccard(distinct_intersect, distinct_union)

        self.assertNotAlmostEqual(similar_jaccard, distinct_jaccard)
