"""Contains a class which describes a network within BiG-SCAPE and methods to help
manipulate this network
"""


# from dependencies
import logging
from typing import Any
from pathlib import Path
from networkx import Graph
from networkx.readwrite import graphml

# from other modules
from src.genbank import BGCRecord
from src.comparison import BGCPair


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

    def write_graphml(self, graph_path: Path):
        """Writes this network graph as a graphml file to the specified output file"""
        graphml.write_graphml_xml(self.graph, graph_path)

    def __contains__(self, __o: Any):
        if not isinstance(__o, BGCPair):
            raise NotImplementedError(
                "Contains check on BSNetwork must be using type BGCPair"
            )

        return self.graph.has_edge(__o.region_a, __o.region_b)
