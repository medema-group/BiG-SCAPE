"""Contains a class which describes a network within BiG-SCAPE and methods to help
manipulate this network
"""


# from dependencies
import logging
from pathlib import Path
from typing import Iterator
from networkx import Graph, connected_components
from networkx.readwrite import graphml
from networkx.readwrite import edgelist

# from other modules
from src.genbank import BGCRecord
from src.comparison import BGCPair

# from this module
from .affinity_propagation import sim_matrix_from_graph, aff_sim_matrix


class BSNetwork:
    """Class to describe a BiG-SCAPE network. Contains methods to add nodes and edges
    with data that is expected in these networks
    """

    def __init__(self):
        self.graph = Graph()

    def add_node(self, node: BGCRecord):
        # the ** before the get_att_dict call converts the dict returned in the
        # method into a keyword argument set
        self.graph.add_node(node, **node.get_attr_dict())

    def add_edge(self, pair: BGCPair, **attr):
        # we want to ensure all edges are accounted for
        if pair.region_a not in self.graph:
            logging.error(
                "Tried to add an edge to the network with a node that does not exist!"
            )
            logging.error("Missing node: %s", pair.region_a)
            raise KeyError(
                "Tried to add an edge to the network with a node that does not exist"
            )

        if pair.region_b not in self.graph:
            logging.error(
                "Tried to add an edge to the network with a node that does not exist!"
            )
            logging.error("Missing node: %s", pair.region_b)
            raise KeyError(
                "Tried to add an edge to the network with a node that does not exist"
            )

        self.graph.add_edge(u_of_edge=pair.region_a, v_of_edge=pair.region_b, **attr)

    def generate_families_cutoff(self, edge_property: str, cutoff: float):
        subgraphs = self.generate_cutoff_subgraphs(edge_property, cutoff)

        family_key = f"family_{cutoff}"

        family_idx = 0
        for subgraph in subgraphs:
            subgraph_dist_matrix = sim_matrix_from_graph(subgraph, edge_property)

            labels, centers = aff_sim_matrix(subgraph_dist_matrix)

            subgraph_family_labels = {}

            for idx, node in enumerate(subgraph.nodes):
                node_family_label = labels[idx]
                if node_family_label not in subgraph_family_labels:
                    subgraph_family_labels[node_family_label] = family_idx
                    family_idx += 1

                node._families[cutoff] = subgraph_family_labels[node_family_label]
                # also add to graph node attributes
                self.graph.nodes.get(node)[family_key] = labels[idx]

    def generate_cutoff_subgraphs(
        self, property_key: str, cutoff: float
    ) -> Iterator[Graph]:
        """Returns an iterator that returns new graphs for each connected component in
        a subgraph created by applying a cutoff for a property on the main graph

        Args:
            property_key (str): the property key to use for the cutoff
            cutoff (float): cutoff value to generate a new subraph for
        """
        # graph.edges is only edge list. no data. use self.graph.adj to get the property
        edge_iter = filter(
            lambda edge: self.graph.adj[edge[0]][edge[1]][property_key] < cutoff,
            self.graph.edges,
        )

        cutoff_subgraph = self.graph.edge_subgraph(edge_iter)

        for connected_component in connected_components(cutoff_subgraph):
            yield self.graph.subgraph(connected_component)

    def write_graphml(self, graph_path: Path):
        """Writes this network graph as a graphml file to the specified output file

        graphml is easier to read in most graph software, but is a lot larger than e.g.
        a TSV file. This makes reading and writing the file more time consuming, and the
        resulting file will occupy more space on your hard drive. Use with care.

        Args:
            graph_path (Path): path to write graphml output to (file)
        """
        graphml.write_graphml_xml(self.graph, graph_path)

    def write_edgelist_tsv(self, graph_path: Path):
        """Writes this network graph as an edge list tsv file to the specified output
        file

        Args:
            graph_path (Path): path to write tsv output to (file)
        """
        fields = ["dist", "jc", "dss", "ai"]
        with open(graph_path, "wb") as tsv_file:
            header = bytes(
                "source\ttarget\t" + "\t".join(fields) + "\n", encoding="utf-8"
            )
            tsv_file.write(header)
            edgelist.write_edgelist(self.graph, tsv_file, data=fields, delimiter="\t")
