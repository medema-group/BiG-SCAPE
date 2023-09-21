"""Contains a class which describes a network within BiG-SCAPE and methods to help
manipulate this network
"""


# from dependencies
import logging
from typing import Any, Optional, Generator
from pathlib import Path
from typing import Iterator
from networkx import Graph, connected_components
from networkx.readwrite import graphml
from networkx.readwrite import edgelist
from sqlalchemy.dialects.sqlite import insert

# from other modules
from src.data import DB
from src.genbank import BGCRecord
from src.comparison import RecordPair
from src.enums import SOURCE_TYPE

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

    def add_edge_pair(self, pair: RecordPair, **attr) -> None:
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

    def get_nodes(
        self, node_types: Optional[list[SOURCE_TYPE]] = None
    ) -> list[BGCRecord]:
        """Returns a list of all nodes in the network.
          If a node type is specified, only node of that kind are returned

        Args:
            node_type (Optional[SOURCE_TYPE], optional): source type of parent gbk. Defaults to None.

        Returns:
            list: nodes

        """
        graph = self.graph
        nodes = list(graph.nodes)

        filtered_nodes = nodes

        if node_types is not None:
            filtered_nodes = []
            for node_type in node_types:
                if not isinstance(node_type, SOURCE_TYPE):
                    raise ValueError("node_type must be of type SOURCE_TYPE")
                nodes_type = filter(
                    lambda node: node.parent_gbk.source_type == node_type, nodes
                )
                filtered_nodes.extend(nodes_type)

        return filtered_nodes

    def get_singletons(
        self, node_types: Optional[list[SOURCE_TYPE]] = None
    ) -> list[BGCRecord]:
        """Returns a list of all nodes in the network that are singletons.

        Args:
            node_types (Optional[list[SOURCE_TYPE]], optional): type of parent gbk. Defaults to None.

        Returns:
            list[BGCRecord]: list of nodes
        """

        singletons = []

        nodes = self.get_nodes(node_types=node_types)

        for node in nodes:
            if list(self.graph.edges(node)) == []:
                singletons.append(node)

        return singletons

    def cull_singletons(self, node_types: Optional[list[SOURCE_TYPE]] = None) -> None:
        """Removes all singletons from the network (of given source type, if specified.)

        Args:
            node_types (Optional[list[SOURCE_TYPE]], optional): parent gbk source type. Defaults to None.
        """
        network = self.graph

        singletons = self.get_singletons(node_types=node_types)

        network.remove_nodes_from(singletons)

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

        for region_a, region_b in self.graph.edges:
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

            DB.execute(upsert_statement, False)

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
        if not isinstance(__o, RecordPair):
            raise NotImplementedError(
                "Contains check on BSNetwork must be using type BGCPair"
            )

        return self.graph.has_edge(__o.region_a, __o.region_b)

    def __eq__(self, __o):
        if not isinstance(__o, BSNetwork):
            return False

        return (
            self.graph.edges == __o.graph.edges and self.graph.nodes == __o.graph.nodes
        )

    @staticmethod
    def load_from_db(bgc_records: set[BGCRecord]):
        """Load a network from the database, using a set of bgc records to determine
        which edges to load
        """
        network = BSNetwork()

        for bgc_record in bgc_records:
            network.add_node(bgc_record)

        distance_table = DB.metadata.tables["distance"]

        record_db_id_dict = {}
        for bgc_record in bgc_records:
            record_db_id_dict[bgc_record._db_id] = bgc_record

        select_statement = distance_table.select().where(
            distance_table.c.region_a_id.in_(record_db_id_dict)
        )

        distance_rows = DB.execute(select_statement).fetchall()

        for distance_row in distance_rows:
            record_a = record_db_id_dict[int(distance_row[0])]
            record_b = record_db_id_dict[int(distance_row[1])]
            network.add_edge(
                record_a,
                record_b,
                dist=distance_row[2],
                jc=distance_row[3],
                ai=distance_row[4],
                dss=distance_row[5],
            )

        return network

    @staticmethod
    def generate_connected_component_networks(
        cutoff: Optional[float] = None,
    ) -> Generator[list[tuple[RecordPair, float, float, float, float]], None, None]:
        """Generate a network for each connected component in the network"""
        # the idea behind this is that the distance table is an edge list. If no
        # distance threshold is applied, any edge is itself a connected component that
        # may be extended with other edges. We can use the distance information that is
        # also present in the network to apply a cutoff value

        if cutoff is None:
            cutoff = 1.0

        # list of nodes to ignore (we have seen them before)
        ignore_nodes = set()

        # get an edge from the database
        edge = get_edge(ignore_nodes)

        # we now have an edge. we need to expand this edge into a connected
        # component. we do this by iteratively selecting for more nodes

        # once edge is none, this means that the get_connected_edges function has
        # found no more new nodes
        while edge is not None:
            # this is our first edge of a connected component
            cc_edges = [edge]

            # we will stop once we have no more new edges
            last_len = len(cc_edges)

            # list of node region ids in the connected component
            edge_nodes = set(cc_edges)

            # we can expand this by adding more edges
            new_edges = get_connected_edges(edge_nodes)
            while len(cc_edges) > last_len:
                edge_nodes.update([edge[0] for edge in cc_edges])
                edge_nodes.update([edge[1] for edge in cc_edges])

                cc_edges.extend(new_edges)

                new_edges = get_connected_edges(edge_nodes)

            # we update the list of ignore nodes with any node we found in the connected
            # component
            ignore_nodes.update([edge[0] for edge in cc_edges])
            ignore_nodes.update([edge[1] for edge in cc_edges])

            # at some point new_edges is empty, so no more new edges were found
            # we can now yield the connected component

            cc = []
            for edge in cc_edges:
                cc.append(RecordPair.load_from_db(edge[0], edge[1]))

            # now we need a new edge to start a new connected component
            edge = get_edge(ignore_nodes)


def get_edge(
    exclude_nodes: set[int],
) -> tuple[RecordPair, float, float, float, float]:
    """Get an edge from the database that is not connected to exclude_nodes"""
    # fetch an edge from the database
    select_statment = (
        DB.metadata.tables["distance"]
        .select()
        .where(DB.metadata.tables["distance"].c.region_a_id.notin_(exclude_nodes))
        .where(DB.metadata.tables["distance"].c.region_b_id.notin_(exclude_nodes))
    )

    edge = DB.execute(select_statment).fetchone()

    return edge


def get_connected_edges(
    include_nodes: set[int],
) -> list[tuple[RecordPair, float, float, float, float]]:
    """Get all edges that are connected to include_nodes, but not to exclude_nodes"""
    # fetch an edge from the database
    select_statment = (
        DB.metadata.tables["distance"]
        .select()
        .where(DB.metadata.tables["distance"].c.region_a_id.in_(include_nodes))
    )

    edges = DB.execute(select_statment).fetchall()

    return edges
