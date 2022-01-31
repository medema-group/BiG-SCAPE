from unittest import TestCase
from src.big_scape.bgc_info import BgcInfo

from src import big_scape

class TestAdjacencyIndex(TestCase):
    """Class containing tests for the run object, which tracks run details and keeps timing
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
