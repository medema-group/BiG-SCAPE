"""Contains functions to manipulate the in-db network"""


# from dependencies
import logging
from typing import Optional, Generator, cast
from sqlalchemy import tuple_, select

# from other modules
from big_scape.data import DB
from big_scape.genbank import BGCRecord


def get_connected_components(
    include_records: list[BGCRecord],
    edge_param_id: int,
    cutoff: Optional[float] = None,
) -> Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
    """Generate a network for each connected component in the network"""
    # the idea behind this is that the distance table is an edge list. If no
    # distance threshold is applied, any edge is itself a connected component that
    # may be extended with other edges. We can use the distance information that is
    # also present in the network to apply a cutoff value

    if cutoff is None:
        cutoff = 1.0

    # list of nodes to include, could be subselection of records in database
    include_nodes: set[int] = set(
        record._db_id for record in include_records if record._db_id is not None
    )

    # list of nodes to ignore (we have seen them before)
    ignore_nodes: set[int] = set()

    # get an edge from the database
    edge = get_edge(include_nodes, ignore_nodes, edge_param_id, cutoff)

    # we now have an edge. we need to expand this edge into a connected
    # component. we do this by iteratively selecting for more nodes

    # once edge is none, this means that the get_connected_edges function has
    # found no more new nodes
    while edge is not None:
        # we can represent our connected component as a set of edges
        connected_component = set((edge,))

        # list of node region ids in the connected component
        edge_node_ids: set[int] = set()
        edge_node_ids.add(edge[0])
        edge_node_ids.add(edge[1])

        # we can expand this by adding more edges
        new_edges = get_connected_edges(
            include_nodes, edge_node_ids, connected_component, edge_param_id, cutoff
        )

        # if we have new edges, we can add them to the connected component
        while len(new_edges) > 0:
            edge_node_ids.update([edge[0] for edge in new_edges])
            edge_node_ids.update([edge[1] for edge in new_edges])

            connected_component.update(new_edges)

            new_edges = get_connected_edges(
                include_nodes, edge_node_ids, connected_component, edge_param_id, cutoff
            )

            # this breaks the while loop if no new edges were found

        # we update the list of ignore nodes with any node we found in the connected
        # component
        ignore_nodes.update([edge[0] for edge in connected_component])
        ignore_nodes.update([edge[1] for edge in connected_component])

        # at some point new_edges is empty, so no more new edges were found
        # we can now yield the connected component

        yield list(connected_component)

        # now we need a new edge to start a new connected component
        edge = get_edge(include_nodes, ignore_nodes, edge_param_id, cutoff)


def get_edge(
    include_nodes: set[int],
    exclude_nodes: set[int],
    edge_param_id: int,
    cutoff: Optional[float] = None,
) -> Optional[tuple[int, int, float, float, float, float, int]]:
    """Get an edge from the database that is not connected to exclude_nodes"""

    if cutoff is None:
        cutoff = 1.0

    # fetch an edge from the database

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]
    select_statment = (
        select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.edge_param_id,
        )
        .where(distance_table.c.record_a_id.in_(include_nodes))
        .where(distance_table.c.record_b_id.in_(include_nodes))
        .where(distance_table.c.record_a_id.notin_(exclude_nodes))
        .where(distance_table.c.record_b_id.notin_(exclude_nodes))
        .where(distance_table.c.edge_param_id == edge_param_id)
        .where(distance_table.c.distance < cutoff)
    )

    edge = DB.execute(select_statment).fetchone()

    if edge is None:
        return None

    return cast(tuple[int, int, float, float, float, float, int], edge)


# TODO: specify edge_param_id?
def get_edges(
    include_nodes: set[int], distance_cutoff: Optional[float] = None
) -> list[tuple[int, int, float, float, float, float, int]]:
    """Get all edges that are connected to include_nodes"""
    if distance_cutoff is None:
        distance_cutoff = 1.0

    # fetch edges from the database

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    distance_table = DB.metadata.tables["distance"]
    select_statement = (
        select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.edge_param_id,
        )
        # equivalent to WHERE (record_a_id in (...) OR record_b_id in (...))
        .where(
            distance_table.c.record_a_id.in_(include_nodes)
            | distance_table.c.record_b_id.in_(include_nodes)
        )
        # equivalent to AND distance < ...
        .where(distance_table.c.distance < distance_cutoff)
    ).compile()

    edges = DB.execute(select_statement).fetchall()

    return cast(list[tuple[int, int, float, float, float, float, int]], edges)


# TODO: check if tested, else test
def get_connected_edges(
    include_nodes: set[int],
    connected_nodes: set[int],
    connected_component: set[tuple[int, int, float, float, float, float, int]],
    edge_param_id: int,
    distance_cutoff: Optional[float] = None,
) -> list[tuple[int, int, float, float, float, float, int]]:
    """Get all edges that are connected to connected_nodes with a certain distance"""
    if distance_cutoff is None:
        distance_cutoff = 1.0

    # fetch edges from the database

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    distance_table = DB.metadata.tables["distance"]
    select_statement = (
        select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.edge_param_id,
        )
        # equivalent to WHERE record_a_id in (...) AND record_b_id in (...)
        .where(distance_table.c.record_a_id.in_(include_nodes))
        .where(distance_table.c.record_b_id.in_(include_nodes))
        # equivalent to AND (record_a_id in (...) OR record_b_id in (...))
        .where(
            distance_table.c.record_a_id.in_(connected_nodes)
            | distance_table.c.record_b_id.in_(connected_nodes)
        )
        # equivalent to AND (record_a_id, record_b_id, ...) NOT IN (connected components)
        .filter(
            ~tuple_(
                distance_table.c.record_a_id,
                distance_table.c.record_b_id,
                distance_table.c.distance,
                distance_table.c.jaccard,
                distance_table.c.adjacency,
                distance_table.c.dss,
                distance_table.c.edge_param_id,
            ).in_(connected_component)
        )
        # equivalent to AND edge_param_id == ...
        .where(distance_table.c.edge_param_id == edge_param_id)
        # equivalent to AND distance < ...
        .where(distance_table.c.distance < distance_cutoff)
        .compile()
    )

    edges = DB.execute(select_statement).fetchall()

    return cast(list[tuple[int, int, float, float, float, float, int]], edges)


def get_query_connected_component(
    include_records: list[BGCRecord],
    query_node_id: Optional[int],
    edge_param_id: int,
    cutoff: Optional[float] = None,
) -> list[tuple[int, int, float, float, float, float, int]]:
    "Generate a network for the query BGC mode connected component in the network"

    if cutoff is None:
        cutoff = 1.0

    # list of nodes to include, could be subselection of records in database
    include_nodes: set[int] = set(
        record._db_id for record in include_records if record._db_id is not None
    )

    # first query edge
    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]

    select_statment = (
        select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.edge_param_id,
        )
        .where(
            (DB.metadata.tables["distance"].c.record_a_id.in_([query_node_id]))
            | (DB.metadata.tables["distance"].c.record_b_id.in_([query_node_id]))
        )
        .where(distance_table.c.edge_param_id == edge_param_id)
        .where(distance_table.c.distance < cutoff)
    )

    edge = DB.execute(select_statment).fetchone()

    if edge is None:
        logging.debug("Query BGC has no edges in database")
        raise RuntimeError("Query BGC has no edges in database")

    # we can represent our connected component as a set of edges
    connected_component = set((edge,))

    # list of node region ids in the connected component
    edge_node_ids: set[int] = set()
    edge_node_ids.add(edge[0])
    edge_node_ids.add(edge[1])

    # we can expand this by adding more edges
    new_edges = get_connected_edges(
        include_nodes, edge_node_ids, connected_component, edge_param_id, cutoff
    )

    # if we have new edges, we can add them to the connected component
    while len(new_edges) > 0:
        edge_node_ids.update([edge[0] for edge in new_edges])
        edge_node_ids.update([edge[1] for edge in new_edges])

        connected_component.update(new_edges)

        new_edges = get_connected_edges(
            include_nodes, edge_node_ids, connected_component, edge_param_id, cutoff
        )

        # this breaks the while loop if no new edges were found

        # at some point new_edges is empty, so no more new edges were found
        # we can now yield the connected component

    return list(connected_component)
