from unittest import TestCase
from src.big_scape.scores import calc_adj_idx
from src.big_scape.bgc_dom_info import BgcDomainInfo
from src.big_scape.distance import calc_ai_pair
from src.big_scape.bgc_info import BgcInfo

from src import big_scape
from test.data.generators import create_cluster_couple

class TestAdjacencyIndex(TestCase):
    """Class containing tests for adjacency index
    """
    test_run: big_scape.Run

    similar_bgc_a: BgcInfo
    similar_bgc_b: BgcInfo

    distinct_bgc_a: BgcInfo
    distinct_bgc_b: BgcInfo

    def setUp(self):
        # similar clusters will have ~80% overlap in domains
        self.similar_bgc_a, self.similar_bgc_b = create_cluster_couple(True, 15)
        
        # for dissimilar we'll do the same but with less overlap
        self.distinct_bgc_a, self.distinct_bgc_b = create_cluster_couple(False, 15)

    def test_adj_index_similar(self):

        cluster_a = self.similar_bgc_a
        cluster_b = self.similar_bgc_b

        a_list = cluster_a.ordered_domain_list
        b_list = cluster_b.ordered_domain_list

        # normally these values are calculated by the BgcDomainInfo class
        # but here we just need to know whether the adjacency index works
        # this is functionally the same as computing the adjacency index for two clusters
        # where there has been no score expansion
        a_start = 0
        b_start = 0
        a_end = len(cluster_a.ordered_domain_list)
        b_end = len(cluster_b.ordered_domain_list)

        ai = calc_adj_idx(a_list, b_list, a_start, a_end, b_start, b_end)

        self.assertNotAlmostEqual(ai, 0.0)

    def test_adj_index_distinct(self):
        cluster_a = self.distinct_bgc_a
        cluster_b = self.distinct_bgc_b

        a_list = cluster_a.ordered_domain_list
        b_list = cluster_b.ordered_domain_list

        # normally these values are calculated by the BgcDomainInfo class
        # but here we just need to know whether the adjacency index works
        # this is functionally the same as computing the adjacency index for two clusters
        # where there has been no score expansion
        a_start = 0
        b_start = 0
        a_end = len(cluster_a.ordered_domain_list)
        b_end = len(cluster_b.ordered_domain_list)

        ai = calc_adj_idx(a_list, b_list, a_start, a_end, b_start, b_end)

        self.assertAlmostEqual(ai, 0.0)

    def test_adj_index_compare(self):
        cluster_a = self.similar_bgc_a
        cluster_b = self.similar_bgc_b

        a_list = cluster_a.ordered_domain_list
        b_list = cluster_b.ordered_domain_list

        a_start = 0
        b_start = 0
        a_end = len(cluster_a.ordered_domain_list)
        b_end = len(cluster_b.ordered_domain_list)

        ai_similar = calc_adj_idx(a_list, b_list, a_start, a_end, b_start, b_end)

        cluster_a = self.distinct_bgc_a
        cluster_b = self.distinct_bgc_b

        a_list = cluster_a.ordered_domain_list
        b_list = cluster_b.ordered_domain_list

        a_start = 0
        b_start = 0
        a_end = len(cluster_a.ordered_domain_list)
        b_end = len(cluster_b.ordered_domain_list)
        
        ai_distinct = calc_adj_idx(a_list, b_list, a_start, a_end, b_start, b_end)

        self.assertNotAlmostEqual(ai_similar, ai_distinct)
