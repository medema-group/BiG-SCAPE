"""Contains utility functions for the comparison module"""

# from python
import logging
import sqlite3

# from dependencies
from sqlalchemy import insert, select

# from other modules
from big_scape.data import DB

# from this module
from .comparable_region import ComparableRegion


def save_edge_to_db(
    edge: tuple[int, int, float, float, float, float, int, ComparableRegion],
    upsert=False,
) -> None:
    """Save edge to the database

    Args:
        edge (tuple[int, int, float, float, float, float, str, int, int, int, int, int,
        int, int, int, bool, ALIGNMENT_MODE,]): edge tuple containing
            record_a_id, record_b_id, distance, jaccard, adjacency, dss, weights,
            lcs start/stop, extension start/stop, reverse, alignment_mode
        upsert (bool, optional): whether to upsert the edge into the database.
    """

    (
        record_a_id,
        record_b_id,
        distance,
        jaccard,
        adjacency,
        dss,
        edge_param_id,
        comparable_region,
    ) = edge

    # save the comparison data to the database

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    distance_table = DB.metadata.tables["distance"]

    # save the entry to the database
    statement = insert(distance_table).values(
        record_a_id=record_a_id,
        record_b_id=record_b_id,
        distance=distance,
        jaccard=jaccard,
        adjacency=adjacency,
        dss=dss,
        edge_param_id=edge_param_id,
        lcs_a_start=comparable_region.lcs_a_start,
        lcs_a_stop=comparable_region.lcs_a_stop,
        lcs_b_start=comparable_region.lcs_b_start,
        lcs_b_stop=comparable_region.lcs_b_stop,
        ext_a_start=comparable_region.a_start,
        ext_a_stop=comparable_region.a_stop,
        ext_b_start=comparable_region.b_start,
        ext_b_stop=comparable_region.b_stop,
        reverse=comparable_region.reverse,
        lcs_domain_a_start=comparable_region.lcs_domain_a_start,
        lcs_domain_a_stop=comparable_region.lcs_domain_a_stop,
        lcs_domain_b_start=comparable_region.lcs_domain_b_start,
        lcs_domain_b_stop=comparable_region.lcs_domain_b_stop,
    )

    if upsert:
        statement = statement.prefix_with("OR REPLACE")

    DB.execute(statement)


def save_edges_to_db(
    edges: list[tuple[int, int, float, float, float, float, ComparableRegion]],
    commit: bool = False,
) -> None:
    """Save many edges to the database

    Args:
        edges (list[tuple[int, int, float, float, float, float, int, int, int, int,
               int, int, int, int, bool]]): list of edges to save
        commit (bool): whether to commit immediately, e.g. during multiprocessing
    """
    # save the comparison data to the database
    # using raw sqlite for this because sqlalchemy is not fast enough

    if not DB.opened():
        raise RuntimeError("DB is not opened")

    if not DB.engine:
        raise RuntimeError("DB.engine is None")

    sqlite_connection = DB.engine.raw_connection().driver_connection

    if not isinstance(sqlite_connection, sqlite3.Connection):
        raise TypeError("Unexpected DB connection type")

    # create a cursor
    cursor = sqlite_connection.cursor()

    # create a query
    # TODO: this should not need ignore. it's there now because protoclusters somehow
    # trigger an integrityerror
    query = "INSERT OR IGNORE INTO distance VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"

    unpacked_edges = [edge[:-1] + edge[-1].to_tuple() for edge in edges]

    cursor.executemany(query, unpacked_edges)

    if commit:
        sqlite_connection.commit()

    # # execute the query. use batches of 100000 to track progress
    # batch_size = 100000
    # with tqdm.tqdm(total=len(edges), unit="edge", desc="Saving edges") as t:
    #     for i in range(0, len(edges), batch_size):
    #         if i + batch_size > len(edges):
    #             batch_size = len(edges) - i

    #         cursor.executemany(query, edges[i : i + batch_size])
    #         sqlite_connection.commit()
    #         t.update(batch_size)


# TODO: does not seem to be used, check if can be removed
# def edges_from_db(
#     pair_generator: RecordPairGenerator,
# ) -> Generator[tuple[RecordPair, float, float, float, float, int], None, None]:
#     """Reconstruct distances from the database instead of recalculating them

#     Args:
#         pair_generator (RecordPairGenerator): pair_generator to regenerate distances for

#     Yields:
#         Generator[tuple[BGCPair, float, float, float, float], None]: generator of distances
#     """
#     # batch size to select database for
#     distance_batch_size = 100

#     # iterate over the pair generator in batches of distance_batch_size
#     for pair_batch in pair_generator.generate_batch(distance_batch_size):
#         # generate a dict of pairs to their ids for record_a and record_b in pair
#         region_ids = {}
#         for pair in pair_batch:
#             region_ids[pair.record_a._db_id] = pair.record_a
#             region_ids[pair.record_b._db_id] = pair.record_b

#         # get the distances from the database

#         if not DB.metadata:
#             raise RuntimeError("DB.metadata is None")

#         distance_table = DB.metadata.tables["distance"]
#         distance_query = distance_table.select().where(
#             distance_table.c.record_a_id.in_(region_ids)
#             & distance_table.c.record_b_id.in_(region_ids)
#         )
#         edges = DB.execute(distance_query).fetchall()

#         # yield the distances
#         for edge in edges:
#             # get the pair
#             record_a = region_ids[edge.record_a_id]
#             record_b = region_ids[edge.record_b_id]
#             pair = RecordPair(record_a, record_b)

#             # get the distances
#             distance: float = edge.distance
#             jaccard: float = edge.jaccard
#             adjacency: float = edge.adjacency
#             dss: float = edge.dss

#             # get the edge param id
#             edge_param_id: int = edge.edge_param_id

#             # yield the distance
#             yield pair, distance, jaccard, adjacency, dss, edge_param_id


def get_edge_param_id(run, weights) -> int:
    """get edge params id if available, else create a new one

    Args:
        run (dict): run parameters
        weights (str): weights category

    Raises:
        RuntimeError: no dabatase

    Returns:
        int: id of the edge param entry
    """

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    alignment_mode = run["alignment_mode"]

    edge_params_table = DB.metadata.tables["edge_params"]
    edge_params_query = (
        select(edge_params_table.c.id)
        .where(edge_params_table.c.alignment_mode == alignment_mode.name)
        .where(edge_params_table.c.weights == weights)
    )

    edge_param_id = DB.execute(edge_params_query).fetchone()

    if edge_param_id is None:
        edge_params_insert = (
            edge_params_table.insert()
            .values(alignment_mode=alignment_mode.name, weights=weights)
            .returning(edge_params_table.c.id)
            .compile()
        )
        cursor_result = DB.execute(edge_params_insert, False)
        edge_param_id = cursor_result.fetchone()

    logging.debug("Edge params id: %d", edge_param_id[0])

    return edge_param_id[0]


def get_edge_weight(edge_param_id: int) -> str:
    """Get edge weights form param ID"""

    if DB.metadata is None:
        raise RuntimeError("DB.metadata is None")

    edge_params_table = DB.metadata.tables["edge_params"]

    edge_weight_query = select(edge_params_table.c.weights).where(
        edge_params_table.c.id == edge_param_id
    )

    weights = DB.execute(edge_weight_query).fetchone()[0]

    return weights
