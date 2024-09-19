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
    delete,
    or_,
    select,
)

# from other modules
from big_scape.data import DB
from big_scape.genbank import BGCRecord
import big_scape.enums as bs_enums
import big_scape.comparison as bs_comparison

# from this module
from big_scape.network.DBAdjList import DBAdjList


def get_connected_components(
    cutoff: float,
    edge_param_id: int,
    bin: bs_comparison.RecordPairGenerator,
    run_id: int,
) -> Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
    """Generate a network for each connected component in the network
        If a seed record is given, the connected component will be generated starting from that record

    Args:
        cutoff (float): the distance cutoff
        edge_param_id (int): the edge parameter id
        bin (bs_comparison.RecordPairGenerator): the bin to generate the connected components for
        run_id (int): current run id

    Yields:
        Generator[list[tuple[int, int, float, float, float, float, int]], None, None]:
        a generator yielding a list of edges for each connected component
    """

    # create a temporary table with the records to include
    include_record_table = None
    if bin is not None:
        include_record_table = create_temp_record_table(bin.source_records)

    # generate connected components using dfs
    generate_connected_components(
        cutoff, edge_param_id, bin.label, run_id, include_record_table
    )

    cc_ids = get_connected_component_ids(
        cutoff,
        bin.label,
        run_id,
        include_record_table,
    )

    logging.info(f"Found {len(cc_ids)} connected components")

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")
    distance_table = DB.metadata.tables["distance"]

    # return connected components per connected component id
    # cc_ids will be repeated accross cutoffs and bins, so
    # we also need to filter by cutoff and bin
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
                        DB.metadata.tables["connected_component"].c.id == cc_id,
                        DB.metadata.tables["connected_component"].c.cutoff == cutoff,
                        DB.metadata.tables["connected_component"].c.bin_label
                        == bin.label,
                        DB.metadata.tables["connected_component"].c.run_id == run_id,
                    )
                ),
                distance_table.c.record_b_id.in_(
                    select(DB.metadata.tables["connected_component"].c.record_id).where(
                        DB.metadata.tables["connected_component"].c.id == cc_id,
                        DB.metadata.tables["connected_component"].c.cutoff == cutoff,
                        DB.metadata.tables["connected_component"].c.bin_label
                        == bin.label,
                        DB.metadata.tables["connected_component"].c.run_id == run_id,
                    )
                ),
                distance_table.c.edge_param_id == edge_param_id,
                distance_table.c.distance < cutoff,
            )
        )

        if include_record_table is not None:
            select_statement = select_statement.where(
                and_(
                    distance_table.c.record_a_id.in_(
                        select(include_record_table.c.record_id)
                    ),
                    distance_table.c.record_b_id.in_(
                        select(include_record_table.c.record_id)
                    ),
                )
            )

        yield list(DB.execute(select_statement))


def dfs(adj_list, start):
    """Depth-first search algorithm

    Args:
        adj_list (dict): adjacency list
        start (int): starting node

    Returns:
        set: set of visited nodes
    """
    stack = [start]
    visited = set()
    while stack:
        node = stack.pop()
        if node not in visited:
            visited.add(node)
            stack.extend([n for n in adj_list[node] if n not in visited])
    return visited


def generate_connected_components(
    cutoff: float,
    edge_param_id: int,
    bin_label: str,
    run_id: int,
    include_record_table: Optional[Table] = None,
    seed_record: Optional[BGCRecord] = None,
) -> None:
    """Generate the connected components for the network with the given parameters

    This uses a rough depth first search to generate the connected components

    If a seed record is given, the connected component will be generated starting from that record
    and the function will return JUST the one connected component

    Args:
        cutoff (Optional[float], optional): the distance cutoff. Defaults to None.
        edge_param_id (int): the edge parameter id
        bin_label (str): the bin label
        run_id (int): the id of the current run
        temp_record_table (Table, optional): a temporary table with the records to include in the
        connected component. Defaults to None.
        seed_record (Optional[BGCRecord], optional): a seed record to start the connected component from.
        Defaults to None.
    """

    db_adj_list = DBAdjList(
        include_record_table,
        cutoff,
        edge_param_id,
        bin_label,
    )

    visited = set()

    if seed_record is not None:
        connected_component = generate_cc_from_node(
            cutoff,
            bin_label,
            run_id,
            db_adj_list,
            seed_record._db_id,
        )
        return

    t = tqdm.tqdm(
        total=len(db_adj_list), unit="nodes", desc="Generating connected components"
    )

    for node in db_adj_list:
        if node not in visited:
            connected_component = generate_cc_from_node(
                cutoff,
                bin_label,
                run_id,
                db_adj_list,
                node,
            )

            visited.update(connected_component)

        t.update(1)

    t.close()


def generate_cc_from_node(
    cutoff: float,
    bin_label: str,
    run_id: int,
    db_adj_list: DBAdjList,
    node: int,
):
    """Generate a connected component from a node using depth first search
    and write it to the database

    Args:
        cutoff (float): the distance cutoff
        bin_label (str): the bin label
        run_id (int): the id of the current run
        db_adj_list (DBAdjList): the adjacency list
        node (int): the starting node

    Returns:
        set: the connected component
    """

    connected_component = dfs(db_adj_list, node)

    cc_id = node

    if len(connected_component) == 1:
        return connected_component

    for node in connected_component:
        DB.execute(
            DB.metadata.tables["connected_component"]
            .insert()
            .values(
                id=cc_id,
                record_id=node,
                cutoff=cutoff,
                bin_label=bin_label,
                run_id=run_id,
            )
        )

    return connected_component


def get_connected_component_ids(
    cutoff: float,
    bin_label: str,
    run_id: int,
    include_record_table: Optional[Table] = None,
) -> list[int]:
    """Get the connected component ids for the given cutoff and edge parameter id

    Args:
        cutoff (float): the distance cutoff
        bin_label (str): label if the current bin
        run_id (int): the id of the current run
        include_record_table (Table, optional): a temporary table with the records to include in the
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
                cc_table.c.bin_label == bin_label,
                cc_table.c.run_id == run_id,
            )
        )
    )

    if include_record_table is not None:
        select_statement = select_statement.where(
            cc_table.c.record_id.in_(select(include_record_table.c.record_id))
        )

    cc_ids = DB.execute(select_statement).fetchall()

    # returned as tuples, convert to list
    cc_ids = [cc_id[0] for cc_id in cc_ids]

    return cc_ids


def get_random_edge(
    cutoff: float,
    edge_param_id: int,
    bin_label: str,
    temp_record_table: Optional[Table] = None,
    seed_record: Optional[BGCRecord] = None,
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

    random_edge_query = (
        # select edge as just record ids
        select(distance_table.c.record_a_id, distance_table.c.record_b_id).where(
            # and where the edge has a distance less than the cutoff and the edge param id is the same
            distance_table.c.distance < cutoff,
            distance_table.c.edge_param_id == edge_param_id,
        )
        # return only one edge
        .limit(1)
    )

    if seed_record is not None:
        random_edge_query = random_edge_query.where(
            or_(
                distance_table.c.record_a_id == seed_record._db_id,
                distance_table.c.record_b_id == seed_record._db_id,
            )
        )

    else:
        random_edge_query = random_edge_query.where(
            # where record a id is not in a connected component with the same cutoff and edge param id
            distance_table.c.record_a_id.notin_(
                select(cc_table.c.record_id).where(
                    cc_table.c.cutoff == cutoff,
                    cc_table.c.bin_label == bin_label,
                )
            )
        ).where(
            # and where record b id is not in a connected component with the same cutoff and edge param id
            distance_table.c.record_b_id.notin_(
                select(cc_table.c.record_id).where(
                    cc_table.c.cutoff == cutoff,
                    cc_table.c.bin_label == bin_label,
                )
            ),
        )

    if temp_record_table is not None:
        random_edge_query = random_edge_query.where(
            and_(
                distance_table.c.record_a_id.in_(select(temp_record_table.c.record_id)),
                distance_table.c.record_b_id.in_(select(temp_record_table.c.record_id)),
            )
        )

    edge = DB.execute(random_edge_query).fetchone()

    return edge


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

    return temp_table


def reset_db_connected_components_table():
    """Removes any data from the connected component table"""
    DB.execute(DB.metadata.tables["connected_component"].delete())


def reference_only_connected_component(connected_component, bgc_records) -> bool:
    """Checks if this connected component is made up of only reference records"""

    has_query = False

    record_ids = []

    for edge in connected_component:
        record_a_id, record_b_id, _, _, _, _, _ = edge
        record_ids.append(record_a_id)
        record_ids.append(record_b_id)

    for record in bgc_records:
        if record._db_id is None:
            raise ValueError("Record has no db id")
        if (
            record._db_id in record_ids
            and record.parent_gbk.source_type == bs_enums.SOURCE_TYPE.QUERY
        ):
            has_query = True
            break

    return not has_query


def get_connected_component_id(
    connected_component: list, cutoff: float, run_id: int
) -> int:
    """Get the connected component id for the given connected component
        expects all edges to be in one connected component, if thats not the
        case, weird things might happen

    Args:
        connected_component: the connected component
        cutoff: the distance cutoff
        run_id: id of the current run

    Returns:
        int: the connected component id
    """

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")

    record_id = connected_component[0][0]

    cc_table = DB.metadata.tables["connected_component"]

    select_statement = (
        select(cc_table.c.id)
        .distinct()
        .where(
            and_(
                cc_table.c.cutoff == cutoff,
                cc_table.c.record_id == record_id,
                cc_table.c.run_id == run_id,
            )
        )
        .limit(1)
    )

    cc_ids = DB.execute(select_statement).fetchone()

    return cc_ids[0]


def remove_connected_component(
    connected_component: list, cutoff: float, run_id: int
) -> None:
    """Removes a connected component from the cc table in the database"""

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")

    cc_id = get_connected_component_id(connected_component, cutoff, run_id)

    cc_table = DB.metadata.tables["connected_component"]

    delete_statement = delete(cc_table).where(
        cc_table.c.id == cc_id, cc_table.c.cutoff == cutoff, cc_table.c.run_id == run_id
    )

    DB.execute(delete_statement)

    DB.commit()
