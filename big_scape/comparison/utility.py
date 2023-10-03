"""Contains utility functions for the comparison module"""

# from python
from typing import Generator

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
