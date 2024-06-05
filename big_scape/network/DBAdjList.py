# from python
from sqlalchemy import Table, select, func
from collections.abc import Mapping

# from bigscape
import big_scape.data as bs_data


class DBAdjList(Mapping):
    """A mapping that represents the adjacency list of a graph stored in the database.

    Args:
        include_records_table (Table): A table that contains the records to include in the graph.
        cutoff (float): The maximum distance between two records for them to be considered adjacent.
        edge_param_id (int): The ID of the edge parameters to use.
        bin_label (str): The label of the bin.
    """

    def __init__(
        self,
        include_records_table: Table,
        cutoff: float,
        edge_param_id: int,
        bin_label: str,
    ):
        self.include_records_table = include_records_table
        self.cutoff = cutoff
        self.edge_param_id = edge_param_id
        self.bin_label = bin_label

        if bs_data.DB.metadata is None:
            bs_data.DB.metadata = bs_data.DB.get_metadata()

        self.bgc_record_table = bs_data.DB.metadata.tables["bgc_record"]
        self.distance_table = bs_data.DB.metadata.tables["distance"]

    def __iter__(self):
        select_query = self.bgc_record_table.select()

        if self.include_records_table is not None:
            select_query = select_query.where(
                self.bgc_record_table.c.id.in_(self.include_records_table)
            )

        cursor_result = bs_data.DB.execute(select_query)

        for record_id in cursor_result.fetchall():
            yield record_id[0]

    def __len__(self):
        if self.include_records_table is None:
            count_query = func.count(self.bgc_record_table.c.id)
        else:
            count_query = func.count(self.include_records_table.c.record_id)

        cursor_result = bs_data.DB.execute(count_query)

        return cursor_result.fetchone()[0]

    def __getitem__(self, key):
        select_a = (
            select(self.distance_table.c.record_a_id)
            .where(self.distance_table.c.record_b_id == key)
            .where(self.distance_table.c.distance < self.cutoff)
            .where(self.distance_table.c.edge_param_id == self.edge_param_id)
        )

        if self.include_records_table is not None:
            select_a = select_a.where(
                self.distance_table.c.record_a_id.in_(self.include_records_table)
            )

        select_b = (
            select(self.distance_table.c.record_b_id)
            .where(self.distance_table.c.record_a_id == key)
            .where(self.distance_table.c.distance < self.cutoff)
            .where(self.distance_table.c.edge_param_id == self.edge_param_id)
        )

        if self.include_records_table is not None:
            select_b = select_b.where(
                self.distance_table.c.record_b_id.in_(self.include_records_table)
            )

        select_query = select_a.union(select_b)

        cursor_result = bs_data.DB.execute(select_query)

        return [record_id[0] for record_id in cursor_result.fetchall()]
