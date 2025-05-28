"""Contains tests for the networking logic of the derep workflow."""

# from python
from unittest import TestCase
from collections import defaultdict
import numpy as np

# from dependencies

# from other modules
from big_scape.dereplicating.networking import Edge, UnionFind, Network


class TestNetworking(TestCase):
    """Test class for the networking logic of the derep workflow."""

    def test_edge_initialization(self):
        """Tests whether the Edge class is initialized correctly."""
        edge = Edge("A", "B", 0.5)
        self.assertEqual(edge.nodeA, "A")
        self.assertEqual(edge.nodeB, "B")
        self.assertEqual(edge.jaccard_similarity, 0.5)

    def test_edge_equality(self):
        """Tests whether two edges are equal regardless of the order of nodes."""
        edge1 = Edge("A", "B", 0.5)
        edge2 = Edge("B", "A", 0.5)
        edge3 = Edge("A", "C", 0.5)

        self.assertEqual(edge1, edge2)
        self.assertNotEqual(edge1, edge3)

    def test_edge_hash(self):
        """Tests whether the hash of an edge is consistent regardless of the order of nodes."""
        edge1 = Edge("A", "B", 0.5)
        edge2 = Edge("B", "A", 0.5)
        edge3 = Edge("A", "C", 0.5)
        edge4 = Edge("A", "B", 0.7)

        self.assertEqual(hash(edge1), hash(edge2))
        self.assertNotEqual(hash(edge1), hash(edge3))
        self.assertNotEqual(hash(edge1), hash(edge4))

    def test_uf_find(self):
        """Tests whether the UnionFind find function correctly finds the representative of a node."""

        # if a node is not yet in the parent dictionary, it
        # will be added with itself as the parent

        uf = UnionFind()
        uf.parent = {"B": "A"}

        self.assertEqual(uf.find("A"), "A")
        self.assertEqual(uf.find("B"), "A")
        self.assertEqual(uf.find("C"), "C")

        expected_parent = {"A": "A", "B": "A", "C": "C"}
        self.assertEqual(uf.parent, expected_parent)

    def test_uf_union(self):
        """Tests whether the UnionFind union function correctly unions two nodes."""

        uf = UnionFind()

        # the second node is processed first, so it will be the parent of the first node
        uf.union("A", "B")

        self.assertEqual(uf.find("A"), "B")
        self.assertEqual(uf.find("B"), "B")

        expected_parent = {"A": "B", "B": "B"}
        self.assertEqual(uf.parent, expected_parent)

        # now C will be the parent of B
        uf.union("B", "C")
        self.assertEqual(uf.find("B"), "C")

        expected_parent = {"A": "B", "B": "C", "C": "C"}
        self.assertEqual(uf.parent, expected_parent)

        # the parent of A is still B since we did not ask for it yet
        # if we ask to find the parent of A, it will get the parent of B which is now C
        # and in the process also set the parent of A to C

        self.assertEqual(uf.find("A"), "C")
        expected_parent = {"B": "C", "A": "C", "C": "C"}
        self.assertEqual(uf.parent, expected_parent)

    def test_network_generate_connected_components(self):
        """Tests whether the generate_connected_components function correctly generates connected components."""

        edges = [
            Edge("A", "B", 0.5),
            Edge("D", "E", 0.5),
            Edge("F", "G", 0.5),
            Edge("X", "X", 1),
        ]

        nodes = set(["A", "B", "D", "E", "F", "G", "X"])

        network = Network(edges, nodes)

        expected_connected_components = {
            "B": {"A", "B"},
            "E": {"D", "E"},
            "G": {"F", "G"},
            "X": {"X"},
        }

        connected_components = network.generate_connected_components()

        self.assertEqual(connected_components, expected_connected_components)

    def test_build_cc_matrix(self):
        """Tests whether the build_cc_matrix function correctly builds a connected component matrix."""

        edges = [
            Edge("A", "B", 0.5),
            Edge("B", "C", 0.7),
            Edge("B", "D", 0.4),
        ]
        nodes = ["A", "B", "C", "D"]

        network = Network(edges, nodes)

        connected_components = network.generate_connected_components()

        # there is only one connected component in this case
        # and the second node on the last edge is the most
        # root parent
        cc_matrix = network.build_cc_matrix(connected_components["D"])

        expected_matrix = np.array(
            [
                # A, B, C, D
                [1.0, 0.5, 0.0, 0.0],  # A
                [0.5, 1.0, 0.7, 0.4],  # B
                [0.0, 0.7, 1.0, 0.0],  # C
                [0.0, 0.4, 0.0, 1.0],  # D
            ]
        )

        self.assertTrue(np.array_equal(cc_matrix, expected_matrix))

    def test_build_cc_matrix_singleton(self):
        """Tests whether the build_cc_matrix function correctly handles a singleton node."""

        edges = [
            Edge("A", "A", 1),
        ]
        nodes = ["A"]

        network = Network(edges, nodes)

        connected_components = network.generate_connected_components()

        cc_matrix = network.build_cc_matrix(connected_components["A"])
        expected_matrix = np.array(
            [
                [1.0],
            ]
        )
        self.assertTrue(np.array_equal(cc_matrix, expected_matrix))

    def test_get_medoid(self):
        """Tests whether the get_medoid function correctly identifies the medoid of a connected component."""

        edges = [
            Edge("A", "B", 0.5),
            Edge("B", "C", 0.7),
            Edge("B", "D", 0.4),
        ]
        nodes = ["A", "B", "C", "D"]
        network = Network(edges, nodes)

        nodes = ["A", "B", "C", "D"]
        expected_matrix = np.array(
            [
                # A, B, C, D
                [1.0, 0.5, 0.0, 0.0],  # A
                [0.5, 1.0, 0.7, 0.4],  # B
                [0.0, 0.7, 1.0, 0.0],  # C
                [0.0, 0.4, 0.0, 1.0],  # D
            ]
        )

        connected_components = network.generate_connected_components()

        medoid = network.get_medoid(connected_components["D"])

        self.assertEqual(medoid, "B")

        expected_index = np.argmax(expected_matrix.sum(axis=0))
        medoid_index = nodes.index(medoid)

        self.assertEqual(expected_index, medoid_index)
        self.assertEqual(medoid, nodes[expected_index])

    def test_get_medoid_singleton(self):
        """Tests whether the get_medoid function correctly identifies the medoid of a singleton connected component."""

        edges_singleton = [Edge("A", "A", 1)]
        nodes_singleton = ["A"]
        network_singleton = Network(edges_singleton, nodes_singleton)
        connected_components_singleton = (
            network_singleton.generate_connected_components()
        )

        medoid_singleton = network_singleton.get_medoid(
            connected_components_singleton["A"]
        )
        self.assertEqual(medoid_singleton, "A")

    def test_set_medoid_representative(self):
        """Tests whether the set_medoid_representative function correctly sets the medoid representative."""

        edges = [
            Edge("A", "B", 0.5),
            Edge("B", "C", 0.7),
            Edge("B", "D", 0.4),
        ]
        nodes = ["A", "B", "C", "D"]
        network = Network(edges, nodes)

        connected_components = network.generate_connected_components()

        rep_connected_components = {}

        # "D" is the parent after the network generation
        network.set_medoid_representative(
            "D", connected_components, rep_connected_components
        )

        expected_representatives = {
            "B": {"A", "B", "C", "D"},  # B is the medoid representative
        }
        self.assertEqual(rep_connected_components, expected_representatives)

    def test_set_medoid_representative_singleton(self):
        """Tests whether the set_medoid_representative function correctly handles a singleton node."""

        edges = [
            Edge("A", "A", 1),
        ]
        nodes = ["A"]

        network = Network(edges, nodes)

        connected_components = network.generate_connected_components()

        rep_connected_components = network.set_medoid_representative(
            "A", connected_components, {}
        )

        expected_representatives = {
            "A": {"A"},
        }
        self.assertEqual(rep_connected_components, expected_representatives)

    def test_build_network(self):
        """Tests whether the build_network function correctly builds the network and generates connected components."""

        # network:
        # A - B   E - F
        #  \ /
        #   C - D  X

        edges = [
            Edge("C", "D", 1),  # CC-A
            Edge("A", "A", 1),  # CC-A
            Edge("A", "B", 1),  # CC-A
            Edge("E", "F", 1),  # CC-B
            Edge("B", "C", 1),  # CC-A
            Edge("X", "X", 1),  # CC-C (singleton)
            Edge("A", "C", 1),  # CC-A
            Edge("B", "B", 1),  # CC-A
        ]
        nodes = ["A", "B", "C", "D", "E", "F", "X"]

        network = Network(edges, nodes)

        rep_connected_components = network.build_network()

        expected_rep_cronnected_components = {
            "C": {"A", "B", "C", "D"},  # CC-A
            "E": {"E", "F"},  # CC-B
            "X": {"X"},  # CC-C (singleton)
        }

        self.assertEqual(rep_connected_components, expected_rep_cronnected_components)

    def test_network_single_cc_linear(self):
        """Tests whether the network correctly identifies a single connected component in a linear graph."""

        # network:
        # A - B - C - D   X - Y  Z

        edges = [
            Edge("C", "D", 1),  # CC-A
            Edge("A", "A", 1),  # CC-A
            Edge("X", "Y", 1),  # CC-B
            Edge("Z", "Z", 1),  # CC-C singleton
            Edge("A", "B", 0.5),  # CC-A
            Edge("B", "C", 1),  # CC-A
            Edge("D", "D", 1),  # CC-A
        ]
        nodes = ["A", "B", "C", "D", "X", "Y", "Z"]
        network = Network(edges, nodes)

        rep_connected_components = network.build_network()

        expected_rep_connected_components = {
            "C": {"A", "B", "C", "D"},  # CC-A
            "X": {"X", "Y"},  # CC-B
            "Z": {"Z"},  # CC-C (singleton)
        }
        self.assertEqual(rep_connected_components, expected_rep_connected_components)
