"""Contains functions to manipulate the in-db network"""


# from dependencies
import logging
from typing import Optional, Generator
from sqlalchemy import tuple_

# from other modules
from src.data import DB
from src.comparison import RecordPair


def get_connected_components(
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
        # we can represent our connected component as a set of edges
        connected_component = set((edge,))

        # list of node region ids in the connected component
        edge_node_ids = set([edge[0], edge[1]])

        # we can expand this by adding more edges
        new_edges = get_connected_edges(edge_node_ids, connected_component, cutoff)

        # if we have new edges, we can add them to the connected component
        while len(new_edges) > 0:
            edge_node_ids.update([edge[0] for edge in connected_component])
            edge_node_ids.update([edge[1] for edge in connected_component])

            connected_component.update(new_edges)

            new_edges = get_connected_edges(edge_node_ids, connected_component, cutoff)

            # this breaks the while loop if no new edges were found

        # we update the list of ignore nodes with any node we found in the connected
        # component
        ignore_nodes.update([edge[0] for edge in connected_component])
        ignore_nodes.update([edge[1] for edge in connected_component])

        # at some point new_edges is empty, so no more new edges were found
        # we can now yield the connected component

        # except if there are only two nodes
        if len(connected_component) == 1:
            edge = get_edge(ignore_nodes)
            continue

        logging.debug(
            "Yielding connected component with %s edges", len(connected_component)
        )

        yield list(connected_component)

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
    connected_component: set[tuple[int, int, float, float, float, float]],
    distance_cutoff: Optional[float] = None,
) -> list[tuple[int, int, float, float, float, float]]:
    """Get all edges that are connected to include_nodes with a certain distance"""
    if distance_cutoff is None:
        distance_cutoff = 1.0

    # fetch an edge from the database
    distance_table = DB.metadata.tables["distance"]
    select_statement = (
        distance_table.select()
        # equivalent to WHERE (region_a_id in (...) OR region_b_id in (...))
        .where(
            distance_table.c.region_a_id.in_(include_nodes)
            | distance_table.c.region_b_id.in_(include_nodes)
        )
        # equivalent to AND (region_a_id, region_b_id, ...) NOT IN (connected components)
        .filter(
            ~tuple_(
                distance_table.c.region_a_id,
                distance_table.c.region_b_id,
                distance_table.c.distance,
                distance_table.c.jaccard,
                distance_table.c.adjacency,
                distance_table.c.dss,
            ).in_(connected_component)
        )
        # equivalent to AND distance < ...
        .where(distance_table.c.distance < distance_cutoff).compile()
    )

    edges = DB.execute(select_statement).fetchall()

    return edges
