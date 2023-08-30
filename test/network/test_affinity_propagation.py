"""Contains tests for affinity propagation"""

# from python
from itertools import combinations
from unittest import TestCase

# from dependencies
from networkx import Graph

# from other modules
from src.network.affinity_propagation import (
    sim_matrix_from_graph,
    aff_sim_matrix,
)


class TestAffinityPropagation(TestCase):
    """Contains tests for affinity propagation"""

    @staticmethod
    def gen_graph():
        """Generates the graph that is used in all tests under this test case"""
        graph = Graph()

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

        for i in range(8):
            graph.add_node(i)

        for a, b in combinations(graph.nodes, 2):
            same_cluster = (a in c1 and b in c1) or (a in c2 and b in c2)
            bridge = sorted([a, b]) in bridges
            if same_cluster:
                distance = 0.01
            elif bridge:
                distance = 0.3
            else:
                # don't add other edges. this is after cutoff
                continue

            graph.add_edge(a, b, dist=distance)

        return graph

    def test_sim_matrix_from_graph(self):
        """Tests whether sim_matrix_from_graph correctly generates a similarity matrix
        from a given graph using the given edge property
        """
        graph = TestAffinityPropagation.gen_graph()

        expected_sim_matrix = [
            [0.0, 0.99, 0.99, 0.99, 0.0, 0.0, 0.0, 0.0],
            [0.99, 0.0, 0.99, 0.99, 0.8, 0.0, 0.0, 0.0],
            [0.99, 0.99, 0.0, 0.99, 0.0, 0.8, 0.0, 0.0],
            [0.99, 0.99, 0.99, 0.0, 0.0, 0.0, 0.8, 0.0],
            [0.0, 0.8, 0.0, 0.0, 0.0, 0.99, 0.99, 0.99],
            [0.0, 0.0, 0.8, 0.0, 0.99, 0.0, 0.99, 0.99],
            [0.0, 0.0, 0.0, 0.8, 0.99, 0.99, 0.0, 0.99],
            [0.0, 0.0, 0.0, 0.0, 0.99, 0.99, 0.99, 0.0],
        ]

        actual_sim_matrix = sim_matrix_from_graph(graph, "dist")

        # convert to list of lists instead of numpy arrays
        actual_sim_matrix = [list(sub_list) for sub_list in expected_sim_matrix]

        self.assertListEqual(expected_sim_matrix, actual_sim_matrix)

    def test_aff_sim_matrix(self):
        """Tests whether affinity propagation is correctly executed on a similarity
        matrix
        """
        graph = TestAffinityPropagation.gen_graph()

        sim_matrix = sim_matrix_from_graph(graph, "dist")

        expected_labels = [0, 0, 0, 0, 1, 1, 1, 1]

        aff_results = aff_sim_matrix(sim_matrix)
        actual_labels = list(aff_results[0])

        self.assertListEqual(expected_labels, actual_labels)
