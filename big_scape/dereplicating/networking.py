"""Module to handle networking of sourmash pairwise distances"""

# from python
from collections import defaultdict
import numpy as np

# from dependencies


# from other modules


class Edge:
    """Class to represent an edge in the network"""

    def __init__(self, nodeA: str, nodeB: str, jaccard_similarity: float) -> None:
        self.nodeA = nodeA
        self.nodeB = nodeB
        self.jaccard_similarity = jaccard_similarity

    def __repr__(self) -> str:
        return f"Edge({self.nodeA}, {self.nodeB}, {self.jaccard_similarity})"

    def __eq__(self, value):
        """Check if two edges are equal, ignoring the order of nodeA and nodeB"""
        # two edges are equal if they have the same nodes, regardless of the order
        # and if they have the same jaccard similarity
        # TODO: is this fine?

        if not isinstance(value, Edge):
            return False
        if self.jaccard_similarity != value.jaccard_similarity:
            return False
        same_nodes = {self.nodeA, self.nodeB} == {value.nodeA, value.nodeB}
        return same_nodes

    def __hash__(self):
        """Hash function for the edge, so it can be used in sets and dicts"""
        # we use a tuple of the nodes and the jaccard similarity as the hash
        # Use frozenset to make the node pair order-insensitive
        return hash((frozenset([self.nodeA, self.nodeB]), self.jaccard_similarity))


class UnionFind:
    """Class to represent a union-find data structure for connected components"""

    def __init__(self):
        self.parent = {}

    def find(self, node):
        """Find the root of the set in which element node is located
        Args:
            node: The element to find the root of
        Returns:
            The root of the set in which element node is located
        """

        # Path compression optimization
        # if node is not its own parent (self.parent.setdefault(x, x)
        # ensures that x is in parent, defaulting to itself.), then recursively
        # find the parent of the parent until we reach the root
        # (i.e. the element that is its own parent)
        if node != self.parent.setdefault(node, node):
            self.parent[node] = self.find(self.parent[node])

        return self.parent[node]

    def union(self, nodeA, nodeB):
        """Union the sets in which elements nodeA and nodeB are located,
        i.e. sets an edge between the two nodes

        Args:
            nodeA: The first element to union
            nodeB: The second element to union
        """

        # Find the roots of the sets in which elements nodeA and nodeB are located
        # and set the parent of the root of nodeA to be the root of nodeB
        # effectively merging the two sets
        self.parent[self.find(nodeA)] = self.find(nodeB)


class Network:
    """Class to represent a network"""

    def __init__(self, edges: list[Edge], nodes: list[str]) -> None:
        self.edges: list[Edge] = edges
        self.nodes: list[str] = nodes
        self.connected_components: dict[str: set[str]] = self.build_network()

    def build_network(self) -> dict[str, set[str]]:
        """Build the network from the edges and generate the connected components

        Returns:
            dict[str, set[str]]: dictionary of connected components, where the key is the medoid representative
            and the value is a set of nodes in the connected component
        """

        # generate the connected components of the network
        connected_components = self.generate_connected_components()

        rep_connected_components = {}
        # for each connected component, we need to set the medoid representative
        # TODO: parallelize this step? what's the best way to do this?
        for parent in connected_components:
            # set the medoid representative for each connected component
            self.set_medoid_representative(
                parent, connected_components, rep_connected_components
            )

        return rep_connected_components

    def generate_connected_components(self) -> defaultdict[str, set[str]]:
        """Generate the connected components of the network using the union-find algorithm

        Returns:
            defaultdict[str, set[str]]: connected component dict, where the key is the parent node,
            and the value is a set of nodes in the connected component
        """

        uf = UnionFind()

        for edge in self.edges:
            uf.union(edge.nodeA, edge.nodeB)

        connected_components = defaultdict(set)
        for node in self.nodes:
            connected_components[uf.find(node)].add(node)

        return connected_components

    def set_medoid_representative(
        self,
        parent: str,
        connected_components: defaultdict[str, set[str]],
        rep_connected_components: dict[str, set[str]],
    ) -> None:
        """Find the medoid of the connected component of the given parent node
        Args:
            parent (str): The parent node to find the medoid for
            connected_components (dict): The connected components of the network
            rep_connected_components (dict): The connected components with medoid representatives
        Returns:
            None: The function modifies the rep_connected_components in place
        """

        if len(connected_components[parent]) == 1:
            # if there is only one node in the connected component, we don't need to do anything
            rep_connected_components[parent] = connected_components[parent]
            return None

        medoid = self.get_medoid(connected_components[parent])

        rep_connected_components[medoid] = connected_components[parent]

        return None

    def get_medoid(self, connected_component: set[str]) -> str:
        """Get the index of the medoid in the similarity matrix

        Args:
            connected_component (set[str]): The connected component to find the medoid for

        Returns:
            str: medoid of the connected component, which is the node with the highest sum of similarities
        """
        # The medoid is the node with the highest sum of similarities

        cc_nodes = sorted(connected_component)
        # get the similarity matrix for the connected component
        similarity_matrix = self.build_cc_matrix(connected_component)

        # get the medoid of the connected component
        medoid_index = np.argmax(similarity_matrix.sum(axis=0))
        medoid = cc_nodes[medoid_index]

        return medoid

    def build_cc_matrix(
        self, connected_component: set[str]
    ) -> np.ndarray:
        """Build a similarity matrix for the connected component of the given parent node
        Args:
            parent (str): The parent node to build the similarity matrix for
        Returns:
            np.ndarray: The similarity matrix for the connected component of the given parent node
        """

        cc_nodes = sorted(connected_component)
        node_to_index = {node: i for i, node in enumerate(cc_nodes)}

        # fill it all with 1.0 so missing edges get similarity of 0
        similarity_matrix = np.full((len(cc_nodes), len(cc_nodes)), 0.0)

        # we're working with similarity so we need to fill the diagonal with 1.0
        np.fill_diagonal(similarity_matrix, 1.0)

        # if there is only one node in the connected component, we don't need to do anything
        if len(cc_nodes) == 1:
            return similarity_matrix

        # the maximum number of edges in a complete graph is n * (n - 1) / 2
        # so we can stop when we reach that number
        max_count_edges = len(cc_nodes) * (len(cc_nodes) - 1) / 2
        edge_count = 0

    # TODO: considerations for performance
    # this will loop through edges for each connected component. this way we can
    # parallelize the building of the similarity matrix for each connected component
    # alternatively, we can loop through the edges just one, and add to the respective
    # matrix. this instead cannot be parallelized.
    # which one is better?
        for edge in self.edges:
            if edge.nodeA in cc_nodes and edge.nodeB in cc_nodes:
                edge_count += 1

                i = node_to_index[edge.nodeA]
                j = node_to_index[edge.nodeB]
                similarity_matrix[i, j] = edge.jaccard_similarity
                similarity_matrix[j, i] = edge.jaccard_similarity

                if edge_count == max_count_edges:
                    break

        return similarity_matrix

    def __repr__(self) -> str:
        """String representation of the network"""
        connected_components = self.connected_components

        # Number of nodes and edges
        num_nodes = len(self.nodes)
        num_edges = len(self.edges)
        num_cc = len(connected_components)

        # Start building the string representation
        repr_str = (
            f"<Network with {num_nodes} nodes, {num_edges} edges, {num_cc} connected components>\n"
            f"Distance matrices for the first 10 components:\n"
        )

        # Add jaccard matrices for a few connected components
        component_count = 0
        for rep in connected_components:
            if component_count >= 10:  # Limit to 10 connected components for brevity
                break
            repr_str += f"Component medoid: {rep}\n"
            nodes = connected_components[rep]
            repr_str += f"Nodes: {', '.join(sorted(nodes))}\n"
            repr_str += "Jaccard similarity matrix:\n"
            matrix = self.build_cc_matrix(nodes)
            repr_str += self.format_matrix(matrix) + "\n"
            component_count += 1

        return repr_str

    def format_matrix(self, matrix: np.ndarray) -> str:
        """
        Format the distance matrix to be a string representation,
        formatting the numbers to 2 decimal places.
        """
        matrix_str = ""
        for row in matrix:
            row_str = "  ".join([f"{x:.1f}" for x in row])
            matrix_str += row_str + "\n"
        return matrix_str
