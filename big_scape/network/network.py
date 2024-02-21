"""Contains functions to manipulate the in-db network"""


# from dependencies
import logging
import tqdm
from typing import Optional, Generator, cast
from sqlalchemy import and_, or_, select, func

# from other modules
from big_scape.data import DB
from big_scape.genbank import BGCRecord


def get_connected_components(
    edge_param_id: int,
    cutoff: Optional[float] = None,
    include_records: Optional[list[BGCRecord]] = None,
) -> Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
    """Generate a network for each connected component in the network"""
    # the idea behind this is that the distance table is an edge list. If no
    # distance threshold is applied, any edge is itself a connected component that
    # may be extended with other edges. We can use the distance information that is
    # also present in the network to apply a cutoff value

    if cutoff is None:
        cutoff = 1.0

    generate_connected_components(edge_param_id, cutoff, include_records)

    cc_ids = get_connected_component_ids(cutoff, edge_param_id, include_records)

    logging.info(f"Found {len(cc_ids)} connected components")

    distance_table = DB.metadata.tables["distance"]

    # return connected components per connected component id
    for cc_id in cc_ids:
        select_statement = select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.edge_param_id,
        ).where(
            and_(
                distance_table.c.record_a_id.in_(
                    select(DB.metadata.tables["connected_component"].c.record_id).where(
                        DB.metadata.tables["connected_component"].c.id == cc_id
                    )
                ),
                distance_table.c.record_b_id.in_(
                    select(DB.metadata.tables["connected_component"].c.record_id).where(
                        DB.metadata.tables["connected_component"].c.id == cc_id
                    )
                ),
                distance_table.c.edge_param_id == edge_param_id,
                distance_table.c.distance < cutoff,
            )
        )

        yield list(DB.execute(select_statement))


def generate_connected_components(
    edge_param_id: int,
    cutoff: float,
    include_records: Optional[list[BGCRecord]] = None,
) -> None:
    """Generate the connected components for the network with the given parameters

    Args:
        include_records (list[BGCRecord]): list of records to include in the connected components
        edge_param_id (int): the edge parameter id
        cutoff (Optional[float], optional): the distance cutoff. Defaults to None.
    """

    distance_table = DB.metadata.tables["distance"]

    edge = get_random_edge(cutoff, edge_param_id, include_records)

    # could be that we already generated all connected components
    if edge is None:
        return

    cc_id = edge[0]

    edges = [edge]

    seen = set()

    cursor = DB.engine.raw_connection().driver_connection.cursor()

    edge_count_query = select(func.count(distance_table.c.record_a_id)).where(
        and_(
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
    )
    num_edges = DB.execute(edge_count_query).fetchone()[0]

    with tqdm.tqdm(total=num_edges) as t:
        while len(edges) > 0:
            inserts = []

            for edge in edges:
                record_id_a = edge[0]
                record_id_b = edge[1]

                if record_id_a not in seen:
                    inserts.append((cc_id, record_id_a, cutoff, edge_param_id))
                if record_id_b not in seen:
                    inserts.append((cc_id, record_id_b, cutoff, edge_param_id))

                seen.add(record_id_a)
                seen.add(record_id_b)

            q = """
            INSERT INTO connected_component (id, record_id, cutoff, edge_param_id)
            VALUES (?, ?, ?, ?);
            """

            cursor.executemany(q, inserts)

            t.update(len(edges))

            last_len = len(edges)

            edges = get_cc_edges(cc_id, cutoff, edge_param_id)

            if len(edges) == last_len:
                # print(f"cc: {', '.join([str(x) for x in seen])}, ({len(seen)})")
                edge = get_random_edge(cutoff, edge_param_id)

                if edge is None:
                    break

                cc_id = edge[0]

                edges = [edge]

                seen = set()

    cursor.close()

    DB.commit()


def get_connected_component_ids(
    cutoff: float, edge_param_id: int, include_records: Optional[list[BGCRecord]] = None
) -> list[int]:
    """Get the connected component ids for the given cutoff and edge parameter id

    Args:
        cutoff (float): the distance cutoff
        edge_param_id (int): the edge parameter id

    Returns:
        list[int]: a list of connected component ids
    """
    cc_table = DB.metadata.tables["connected_component"]
    select_statement = (
        select(cc_table.c.id)
        .distinct()
        .where(
            and_(
                cc_table.c.cutoff == cutoff,
                cc_table.c.edge_param_id == edge_param_id,
            )
        )
    )

    if include_records is not None:
        include_ids = [record._db_id for record in include_records]
        select_statement = select_statement.where(cc_table.c.record_id.in_(include_ids))

    cc_ids = DB.execute(select_statement).fetchall()

    # returned as tuples, convert to list
    cc_ids = [cc_id[0] for cc_id in cc_ids]

    return cc_ids


def get_random_edge(
    cutoff: float, edge_param_id: int, include_records: Optional[list[BGCRecord]] = None
) -> Optional[tuple[int, int]]:
    """
    Get a random edge from the database that is not in any connected component
    and has a distance less than the cutoff

    Note that this returns only the ids to reduce the amount of data

    Args:
        include_records: the records to include in the connected component
        cutoff: the distance cutoff
        edge_param_id: the edge parameter id

    Returns:
        Optional[tuple[int, int]]: a tuple with the record ids of the edge or none
    """

    distance_table = DB.metadata.tables["distance"]
    cc_table = DB.metadata.tables["connected_component"]

    random_edge_query = (
        select(distance_table.c.record_a_id, distance_table.c.record_b_id)
        .where(
            and_(
                distance_table.c.record_a_id.notin_(select(cc_table.c.record_id)),
                distance_table.c.record_b_id.notin_(select(cc_table.c.record_id)),
                distance_table.c.distance < cutoff,
                distance_table.c.edge_param_id == edge_param_id,
            )
        )
        .limit(1)
    )

    if include_records is not None:
        include_ids = [record._db_id for record in include_records]
        random_edge_query = random_edge_query.where(
            or_(
                distance_table.c.record_a_id.in_(include_ids),
                distance_table.c.record_b_id.in_(include_ids),
            )
        )

    edge = DB.execute(random_edge_query).fetchone()

    return edge


def get_cc_edges(
    cc_id: int, cutoff: float, edge_param_id: int
) -> list[tuple[int, int]]:
    """
    Get all edges connected to an existing connected component with a distance less than
    the cutoff

    Args:
        cc_id: the connected component id
        cutoff: the distance cutoff

    Returns:
        Optional[tuple[int, int]]: a tuple with the record ids of the edge or none
    """

    distance_table = DB.metadata.tables["distance"]
    cc_table = DB.metadata.tables["connected_component"]

    cc_edge_query = select(
        distance_table.c.record_a_id, distance_table.c.record_b_id
    ).where(
        and_(
            or_(
                distance_table.c.record_a_id.in_(
                    select(cc_table.c.record_id).where(cc_table.c.id == cc_id)
                ),
                distance_table.c.record_b_id.in_(
                    select(cc_table.c.record_id).where(cc_table.c.id == cc_id)
                ),
            ),
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
    )
    edges = DB.execute(cc_edge_query).fetchall()

    return edges


# generic functions


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


def get_nodes_from_cc(
    connected_component, bgc_records: list[BGCRecord]
) -> list[BGCRecord]:
    "get the nodes from the connected component"

    cc_record_ids = set()
    cc_record_list = []

    for edge in connected_component:
        (
            record_a_id,
            record_b_id,
            dist,
            jacc,
            adj,
            dss,
            edge_param_id,
        ) = edge
        cc_record_ids.add(record_a_id)
        cc_record_ids.add(record_b_id)

    for record in bgc_records:
        if record._db_id is None:
            raise ValueError("Region in bin has no db id!")
        if record._db_id not in cc_record_ids:
            continue

        cc_record_list.append(record)

    return cc_record_list
