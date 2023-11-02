"""Contains tests for affinity propagation"""

# from python
from itertools import combinations
from unittest import TestCase

# from other modules
import big_scape.network.families as bs_families


class TestAffinityPropagation(TestCase):
    """Contains tests for affinity propagation"""

    @staticmethod
    def gen_edge_list():
        """Generates the edge list that is used in all tests under this test case"""

        # graph topology:
        # (numbers indicate distance)
        #    1---0.3---4
        #   /|\       /|\
        #  0-+-2-0.3-5-+-7
        #   \|/       \|/
        #    3---0.3---6
        #  |---|     |---|
        # cluster1  cluster2
        # within-cluster distances are set to 0.01

        c1 = [0, 1, 2, 3]
        c2 = [4, 5, 6, 7]
        bridges = [
            [1, 4],
            [2, 5],
            [3, 6],
        ]

        edge_list = []

        for a, b in combinations(range(8), 2):
            same_cluster = (a in c1 and b in c1) or (a in c2 and b in c2)
            bridge = sorted([a, b]) in bridges
            if same_cluster:
                distance = 0.01
            elif bridge:
                distance = 0.3
            else:
                # don't add other edges. this is after cutoff
                continue

            edge_list.append((a, b, distance, 0, 0, 0, ""))

        return edge_list

    def test_sim_matrix_from_graph(self):
        """Tests whether sim_matrix_from_graph correctly generates a similarity matrix
        from a given graph using the given edge property
        """
        adj_list = TestAffinityPropagation.gen_edge_list()

        expected_sim_matrix = [
            [1.0, 0.99, 0.99, 0.99, 0.0, 0.0, 0.0, 0.0],
            [0.99, 1.0, 0.99, 0.99, 0.7, 0.0, 0.0, 0.0],
            [0.99, 0.99, 1.0, 0.99, 0.0, 0.7, 0.0, 0.0],
            [0.99, 0.99, 0.99, 1.0, 0.0, 0.0, 0.7, 0.0],
            [0.0, 0.7, 0.0, 0.0, 1.0, 0.99, 0.99, 0.99],
            [0.0, 0.0, 0.7, 0.0, 0.99, 1.0, 0.99, 0.99],
            [0.0, 0.0, 0.0, 0.7, 0.99, 0.99, 1.0, 0.99],
            [0.0, 0.0, 0.0, 0.0, 0.99, 0.99, 0.99, 1.0],
        ]

        actual_sim_matrix, _ = bs_families.edge_list_to_sim_matrix(adj_list)

        actual_sim_matrix = actual_sim_matrix.tolist()

        # check if the matrices have the exact same values
        self.assertListEqual(expected_sim_matrix, actual_sim_matrix)

    def test_aff_sim_matrix(self):
        """Tests whether affinity propagation is correctly executed on a similarity
        matrix
        """
        adj_list = TestAffinityPropagation.gen_edge_list()

        sim_matrix, _ = bs_families.edge_list_to_sim_matrix(adj_list)

        expected_labels = [0, 0, 0, 0, 1, 1, 1, 1]

        aff_results = bs_families.aff_sim_matrix(sim_matrix)
        actual_labels = list(aff_results[0])

        self.assertListEqual(expected_labels, actual_labels)
