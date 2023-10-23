"""Contains utility functions for the comparison module"""

# from python
from typing import Generator, cast
import sqlite3

# from dependencies
from sqlalchemy import insert

# from other modules
from big_scape.data import DB
from big_scape.comparison.binning import RecordPairGenerator, RecordPair


def save_edge_to_db(
    edge: tuple[int, int, float, float, float, float], upsert=False
) -> None:
    """Save edge to the database

    Args:
        edge (tuple[int, int, float, float, float, float]): edge tuple containing
            region_a_id, region_b_id, distance, jaccard, adjacency, dss
        upsert (bool, optional): whether to upsert the edge into the database.
    """

    region_a_id, region_b_id, distance, jaccard, adjacency, dss = edge

    # save the comparison data to the database

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    distance_table = DB.metadata.tables["distance"]

    # save the entry to the database
    statement = insert(distance_table).values(
        region_a_id=region_a_id,
        region_b_id=region_b_id,
        distance=distance,
        jaccard=jaccard,
        adjacency=adjacency,
        dss=dss,
    )

    if upsert:
        statement = statement.prefix_with("OR REPLACE")

    DB.execute(statement)


def save_edges_to_db(edges: list[tuple[int, int, float, float, float, float]]) -> None:
    """Save many edges to the database

    Args:
        edges (list[tuple[int, int, float, float, float, float]]): list of edges to save
    """
    # save the comparison data to the database
    # using raw sqlite for this because sqlalchemy is not fast enough

    if not DB.opened():
        raise RuntimeError("DB is not opened")

    if not DB.engine:
        raise RuntimeError("DB.engine is None")

    sqlite_connection = DB.engine.raw_connection()
    sqlite_connection = cast(sqlite3.Connection, sqlite_connection)

    # create a cursor
    cursor = sqlite_connection.cursor()

    # create a query
    # TODO: this should not need ignore. it's there now because protoclusters somehow
    # trigger an integrityerror
    query = "INSERT OR IGNORE INTO distance VALUES (?, ?, ?, ?, ?, ?)"

    cursor.executemany(query, edges)

    # # execute the query. use batches of 100000 to track progress
    # batch_size = 100000
    # with tqdm.tqdm(total=len(edges), unit="edge", desc="Saving edges") as t:
    #     for i in range(0, len(edges), batch_size):
    #         if i + batch_size > len(edges):
    #             batch_size = len(edges) - i

    #         cursor.executemany(query, edges[i : i + batch_size])
    #         sqlite_connection.commit()
    #         t.update(batch_size)


def edges_from_db(
    pair_generator: RecordPairGenerator,
) -> Generator[tuple[RecordPair, float, float, float, float], None, None]:
    """Reconstruct distances from the database instead of recalculating them

    Args:
        pair_generator (RecordPairGenerator): pair_generator to regenerate distances for

    Yields:
        Generator[tuple[BGCPair, float, float, float, float]]: generator of distances
    """
    # batch size to select database for
    distance_batch_size = 100

    # iterate over the pair generator in batches of distance_batch_size
    for pair_batch in pair_generator.generate_batch(distance_batch_size):
        # generate a dict of pairs to their ids for region_a and region_b in pair
        region_ids = {}
        for pair in pair_batch:
            region_ids[pair.region_a._db_id] = pair.region_a
            region_ids[pair.region_b._db_id] = pair.region_b

        # get the distances from the database

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]
        distance_query = distance_table.select().where(
            distance_table.c.region_a_id.in_(region_ids)
            & distance_table.c.region_b_id.in_(region_ids)
        )
        edges = DB.execute(distance_query).fetchall()

        # yield the distances
        for edge in edges:
            # get the pair
            region_a = region_ids[edge.region_a_id]
            region_b = region_ids[edge.region_b_id]
            pair = RecordPair(region_a, region_b)

            # get the distances
            distance: float = edge.distance
            jaccard: float = edge.jaccard
            adjacency: float = edge.adjacency
            dss: float = edge.dss

            # yield the distance
            yield pair, distance, jaccard, adjacency, dss
