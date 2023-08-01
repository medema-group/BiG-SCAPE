"""Contains a class which describes a network within BiG-SCAPE and methods to help
manipulate this network
"""


# from dependencies
import logging
from typing import Any
from pathlib import Path
from typing import Iterator
from networkx import Graph, connected_components
from networkx.readwrite import graphml
from networkx.readwrite import edgelist
from sqlalchemy.dialects.sqlite import insert

# from other modules
from src.data import DB
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

    def add_node(self, node: BGCRecord) -> None:
        """Add a node in the form of a BGCRecord. This function calls node.get_attr_dict() to get
        the attributes that are relevant for this node (e.g. family memberships)

        Args:
            node (BGCRecord): BGCRecord
        """
        # the ** before the get_att_dict call converts the dict returned in the
        # method into a keyword argument set
        self.graph.add_node(node, **node.get_attr_dict())

    def add_edge_pair(self, pair: BGCPair, **attr) -> None:
        """Add an edge between two regions in a pair.

        Args:
            pair (BGCPair): BGC pair to add an edge for
            attr (dict[str, Any]): values to add to edge
        """
        self.add_edge(pair.region_a, pair.region_b, **attr)

    def add_edge(self, region_a: BGCRecord, region_b: BGCRecord, **attr) -> None:
        """Add an edge between two regions

        Args:
            region_a (BGCRecord): U node of edge as a BGCRecord
            region_b (BGCRecord): V node of edge as a BGCRecord

        Raises:
            KeyError: Raised if a node from an edge is not present in the network
        """
        # we want to ensure all edges are accounted for
        if region_a not in self.graph:
            logging.error(
                "Tried to add an edge to the network with a node that does not exist!"
            )
            logging.error("Missing node: %s", region_a)
            raise KeyError(
                "Tried to add an edge to the network with a node that does not exist"
            )

        if region_b not in self.graph:
            logging.error(
                "Tried to add an edge to the network with a node that does not exist!"
            )
            logging.error("Missing node: %s", region_b)
            raise KeyError(
                "Tried to add an edge to the network with a node that does not exist"
            )

        self.graph.add_edge(u_of_edge=region_a, v_of_edge=region_b, **attr)

    def generate_families_cutoff(self, edge_property: str, cutoff: float) -> None:
        """Generate the families for nodes in a network, using a given cutoff for a property that is
        present on an edge. This will usually be 'dist'

        Args:
            edge_property (str): edge property to use for cutoff
            cutoff (float): the cutoff value to use
        """
        subgraphs = self.generate_cutoff_subgraphs(edge_property, cutoff)

        family_key = f"family_{cutoff}"

        # init family centers
        for idx, node in enumerate(self.graph.nodes):
            node._families[cutoff] = idx

        for subgraph in subgraphs:
            # do not cluster small subgraphs
            if len(subgraph) < 3:
                for node in subgraph.nodes:
                    node._families[cutoff] = list(subgraph.nodes)[0]._families[cutoff]
                continue

            subgraph_dist_matrix = sim_matrix_from_graph(subgraph, edge_property)

            labels, centers = aff_sim_matrix(subgraph_dist_matrix)

            for idx, node in enumerate(subgraph.nodes):
                node._families[cutoff] = centers[labels[idx]]
                # also add to graph node attributes
                self.graph.nodes.get(node)[family_key] = centers[labels[idx]]

    def export_distances_to_db(self, commit=True) -> None:
        """Save the distances stored in this network to the database.

        If a distance between two nodes already exists, this will overwrite that
        distance

        Args:
            commit (bool, optional): Commit the result immediately. Defaults to True.
        """
        distance_table = DB.metadata.tables["distance"]

        for region_a in self.graph.adj:
            for region_b in self.graph.adj[region_a]:
                edge_data = self.graph.adj[region_a][region_b]

                upsert_statement = (
                    insert(distance_table)
                    .values(
                        region_a_id=region_a._db_id,
                        region_b_id=region_b._db_id,
                        distance=edge_data["dist"],
                        jaccard=edge_data["jc"],
                        adjacency=edge_data["ai"],
                        dss=edge_data["dss"],
                    )
                    .on_conflict_do_update(
                        index_elements=["region_a_id", "region_b_id"],
                        set_={
                            "distance": edge_data["dist"],
                            "jaccard": edge_data["jc"],
                            "adjacency": edge_data["ai"],
                            "dss": edge_data["dss"],
                        },
                    )
                )

                DB.execute(upsert_statement)

        if commit:
            DB.commit()

    def generate_cutoff_subgraphs(
        self, property_key: str, cutoff: float
    ) -> Iterator[Graph]:
        """Return an iterator that returns new graphs for each connected component in
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

        cutoff_subgraph = Graph()
        cutoff_subgraph.add_nodes_from(self.graph)

        cutoff_subgraph.add_edges_from(edge_iter)

        # cutoff_subgraph = self.graph.edge_subgraph(edge_iter)

        for connected_component in connected_components(cutoff_subgraph):
            yield self.graph.subgraph(connected_component)

    def write_graphml(self, graph_path: Path) -> None:
        """Write this network graph as a graphml file to the specified output file

        graphml is easier to read in most graph software, but is a lot larger than e.g.
        a TSV file. This makes reading and writing the file more time consuming, and the
        resulting file will occupy more space on your hard drive. Use with care.

        Args:
            graph_path (Path): path to write graphml output to (file)
        """
        graphml.write_graphml_xml(self.graph, graph_path)

    def write_edgelist_tsv(self, graph_path: Path) -> None:
        """Write this network graph as an edge list tsv file to the specified output
        file

        Args:
            graph_path (Path): path to write tsv output to (file)
        """
        fields = ["dist", "jc", "ai", "dss"]
        with open(graph_path, "wb") as tsv_file:
            header = bytes(
                "source\ttarget\t" + "\t".join(fields) + "\n", encoding="utf-8"
            )
            tsv_file.write(header)
            edgelist.write_edgelist(self.graph, tsv_file, data=fields, delimiter="\t")

    def __contains__(self, __o: Any):
        if not isinstance(__o, BGCPair):
            raise NotImplementedError(
                "Contains check on BSNetwork must be using type BGCPair"
            )

        return self.graph.has_edge(__o.region_a, __o.region_b)
