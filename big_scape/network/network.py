"""Contains functions to manipulate the in-db network"""


# from dependencies
import logging
import random
import string
import tqdm
from typing import Optional, Generator, cast
from sqlalchemy import (
    Column,
    ForeignKey,
    Integer,
    Table,
    and_,
    distinct,
    or_,
    select,
    func,
)

# from other modules
from big_scape.data import DB
from big_scape.genbank import BGCRecord


def get_connected_components(
    cutoff: float,
    edge_param_id: int,
    include_records: Optional[list[BGCRecord]] = None,
) -> Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
    """Generate a network for each connected component in the network

    Args:
        cutoff (float): the distance cutoff
        edge_param_id (int): the edge parameter id
        include_records (list[BGCRecord]], Optional): list of records to include in
        the connected components. Defaults to None.

    Yields:
        Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
        a generator yielding a list of edges for each connected component
    """

    # create a temporary table with the records to include
    temp_record_table = None
    if include_records is not None:
        temp_record_table = create_temp_record_table(include_records)

    # generate connected components using dfs
    generate_connected_components(cutoff, edge_param_id, temp_record_table)

    cc_ids = get_connected_component_ids(cutoff, edge_param_id, temp_record_table)

    logging.info(f"Found {len(cc_ids)} connected components")

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
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
    cutoff: float, edge_param_id: int, temp_record_table: Optional[Table] = None
) -> None:
    """Generate the connected components for the network with the given parameters

    This uses a rough depth first search to generate the connected components

    Args:
        cutoff (Optional[float], optional): the distance cutoff. Defaults to None.
        edge_param_id (int): the edge parameter id
        temp_record_table (Table, optional): a temporary table with the records to include in the
        connected component. Defaults to None.
    """

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]

    edge = get_random_edge(cutoff, edge_param_id, temp_record_table)

    # could be that we already generated all connected components
    if edge is None:
        return

    cc_id = edge[0]

    edges = [edge]

    seen = set()

    if DB.engine is None:
        raise RuntimeError("DB.engine is None")
    cursor = DB.engine.raw_connection().driver_connection.cursor()

    edge_count_query = select(func.count(distance_table.c.record_a_id)).where(
        and_(
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
    )
    num_edges = DB.execute(edge_count_query).fetchone()[0]

    with tqdm.tqdm(total=num_edges, desc="Generating connected components") as t:
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

            last_len = len(edges)

            edges = get_cc_edges(cc_id, cutoff, edge_param_id)

            if last_len == len(edges):
                # print(f"cc: {', '.join([str(x) for x in seen])}, ({len(seen)})")
                edge = get_random_edge(cutoff, edge_param_id, temp_record_table)

                t.update(len(edges))

                if edge is None:
                    break

                cc_id = edge[0]

                edges = [edge]

                seen = set()

    cursor.close()

    DB.commit()


# TODO: not used, delete
def has_missing_cc_assignments(
    cutoff: float, edge_param_id: int, temp_record_table: Optional[Table] = None
) -> bool:
    """Check if there are any missing connected component assignments for the given cutoff and edge parameter id

    Args:
        cutoff (float): the distance cutoff
        edge_param_id (int): the edge parameter id
        temp_table (Table, optional): a temporary table with the records to include in the connected
        component. Defaults to None.

    Returns:
        bool: True if there are missing connected component assignments, False otherwise
    """

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]
    cc_table = DB.metadata.tables["connected_component"]

    select_statement = (
        select(func.count(distinct(distance_table.c.record_a_id)))
        .where(
            and_(
                distance_table.c.distance < cutoff,
                distance_table.c.edge_param_id == edge_param_id,
            )
        )
        .where(distance_table.c.record_a_id.notin_(select(cc_table.c.record_id)))
    )

    if temp_record_table is not None:
        select_statement = select_statement.where(
            distance_table.c.record_a_id.in_(select(temp_record_table.c.record_id))
        )

    num_missing = DB.execute(select_statement).fetchone()[0]

    return num_missing > 0


def get_connected_component_ids(
    cutoff: float, edge_param_id: int, temp_record_table: Optional[Table] = None
) -> list[int]:
    """Get the connected component ids for the given cutoff and edge parameter id

    Args:
        cutoff (float): the distance cutoff
        edge_param_id (int): the edge parameter id
        temp_record_table (Table, optional): a temporary table with the records to include in the
        connected component. Defaults to None.

    Returns:
        list[int]: a list of connected component ids
    """
    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
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

    if temp_record_table is not None:
        select_statement = select_statement.where(
            cc_table.c.record_id.in_(select(temp_record_table.c.record_id))
        )

    cc_ids = DB.execute(select_statement).fetchall()

    # returned as tuples, convert to list
    cc_ids = [cc_id[0] for cc_id in cc_ids]

    return cc_ids


def get_random_edge(
    cutoff: float, edge_param_id: int, temp_record_table: Optional[Table] = None
) -> Optional[tuple[int, int]]:
    """
    Get a random edge from the database that is not in any connected component
    and has a distance less than the cutoff

    Note that this returns only the ids to reduce the amount of data

    Args:
        cutoff: the distance cutoff
        edge_param_id: the edge parameter id
        temp_record_table (Table, optional): a temporary table with the records to include in the
        connected component. Defaults to None.

    Returns:
        Optional[tuple[int, int]]: a tuple with the record ids of the edge or None
    """
    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]
    cc_table = DB.metadata.tables["connected_component"]

    # this query is complicated, breaking it down:

    random_edge_query = (
        # select edge as just record ids
        select(distance_table.c.record_a_id, distance_table.c.record_b_id)
        .where(
            # where record a id is not in a connected component with the same cutoff and edge param id
            distance_table.c.record_a_id.notin_(
                select(cc_table.c.record_id).where(
                    and_(
                        cc_table.c.cutoff == cutoff,
                        cc_table.c.edge_param_id == edge_param_id,
                    )
                )
            )
        )
        .where(
            # and where record b id is not in a connected component with the same cutoff and edge param id
            distance_table.c.record_b_id.notin_(
                select(cc_table.c.record_id).where(
                    and_(
                        cc_table.c.cutoff == cutoff,
                        cc_table.c.edge_param_id == edge_param_id,
                    )
                )
            ),
        )
        .where(
            # and where the edge has a distance less than the cutoff and the edge param id is the same
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
        # return only one edge
        .limit(1)
    )

    if temp_record_table is not None:
        random_edge_query = random_edge_query.where(
            or_(
                distance_table.c.record_a_id.in_(select(temp_record_table.c.record_id)),
                distance_table.c.record_b_id.in_(select(temp_record_table.c.record_id)),
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
    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]
    cc_table = DB.metadata.tables["connected_component"]

    cc_edge_query = select(
        distance_table.c.record_a_id, distance_table.c.record_b_id
    ).where(
        and_(
            or_(
                distance_table.c.record_a_id.in_(
                    select(cc_table.c.record_id)
                    .where(cc_table.c.id == cc_id)
                    .where(cc_table.c.cutoff == cutoff)
                ),
                distance_table.c.record_b_id.in_(
                    select(cc_table.c.record_id)
                    .where(cc_table.c.id == cc_id)
                    .where(cc_table.c.cutoff == cutoff)
                ),
            ),
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
    )
    edges = DB.execute(cc_edge_query).fetchall()

    return edges


# generic functions
# TODO: not used, delete
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


def create_temp_record_table(include_records: list[BGCRecord]) -> Table:
    """Create a temporary table with ids of given records

    Args:
        include_records (list[BGCRecord]): the records to include in the connected component

    Returns:
        Table: the temporary table
    """
    ids = []

    for record in include_records:
        if record._db_id is None:
            raise ValueError("Record has no db id")
        ids.append(record._db_id)

    # generate a short random string
    temp_table_name = "temp_" + "".join(random.choices(string.ascii_lowercase, k=10))

    temp_table = Table(
        temp_table_name,
        DB.metadata,
        Column(
            "record_id",
            Integer,
            ForeignKey(DB.metadata.tables["bgc_record"].c.id),
            primary_key=True,
            nullable=False,
        ),
        prefixes=["TEMPORARY"],
    )

    DB.metadata.create_all(DB.engine)

    # create_temp_table = f"""
    #     CREATE TEMPORARY TABLE {temp_table_name} (
    #         record_id INTEGER PRIMARY KEY NOT NULL references bgc_record(id)
    #     );
    # """

    # DB.execute_raw_query(create_temp_table)

    if DB.engine is None:
        raise RuntimeError("DB engine is None")

    cursor = DB.engine.raw_connection().driver_connection.cursor()

    insert_query = f"""
        INSERT INTO {temp_table_name} (record_id) VALUES (?);
    """

    cursor.executemany(insert_query, [(x,) for x in ids])  # type: ignore

    cursor.close()

    DB.commit()

    if DB.metadata is None:
        raise ValueError("DB metadata is None")

    table = Table(temp_table_name, DB.metadata, autoload_with=DB.engine)

    return table


def reset_db_connected_components():
    """Removes any data from the connected component table"""
    DB.execute(DB.metadata.tables["connected_component"].delete())
