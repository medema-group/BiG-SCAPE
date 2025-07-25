"""Contains classes and functions for generating bins of Regions to compare

At this level, the comparisons are referred to as pairs. Whenever anything talks about
pairs, it refers to things generated from these classes. This is distinct from what are
referred to as edges, which are pairs that have a (set of) distances between them and
may be present in the database.

TODO: this file is very long and is begging for refactoring. a lot of the classes
are very similar apart from the method in which they query the database. This could
almost certainly be abstracted somehow
"""

# from python
from __future__ import annotations
import logging
from itertools import combinations
import random
import string
from typing import Generator, Iterator, Optional
from sqlalchemy import Column, ForeignKey, Integer, Table, and_, select, func, or_

# from other modules
from big_scape.cli.config import BigscapeConfig
from big_scape.data import DB
from big_scape.genbank import (
    BGCRecord,
    Region,
    ProtoCluster,
    ProtoCore,
)
from big_scape.enums import SOURCE_TYPE, CLASSIFY_MODE, RECORD_TYPE

import big_scape.comparison as bs_comparison


def create_temp_record_id_table(record_ids: set[int]) -> Table:
    """Create a temporary table with ids of given records

    Args:
        record_ids (list[int]): the ids of the records to add to the temporary table

    Returns:
        Table: the temporary table
    """

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

    if DB.engine is None:
        raise RuntimeError("DB engine is None")

    cursor = DB.engine.raw_connection().driver_connection.cursor()

    insert_query = f"""
        INSERT INTO {temp_table_name} (record_id) VALUES (?);
    """

    # local function for batching
    def batch_hash(record_ids: list[int], n: int):
        count = 0
        batch = []
        for i in record_ids:
            batch.append(i)
            count += 1
            if count % n == 0:
                yield batch
                batch = []

        if len(batch) > 0:
            yield batch

    for hash_batch in batch_hash(record_ids, 1000):
        cursor.executemany(insert_query, [(x,) for x in hash_batch])  # type: ignore

    cursor.close()

    DB.commit()

    if DB.metadata is None:
        raise ValueError("DB metadata is None")

    return temp_table


# weights are in the order JC, AI, DSS, Anchor boost
LEGACY_WEIGHTS = {
    "PKSI": {"weights": (0.22, 0.02, 0.76, 1.0)},
    "PKSother": {"weights": (0.0, 0.68, 0.32, 4.0)},
    "NRPS": {"weights": (0.0, 0.0, 1.0, 4.0)},
    "RiPP": {"weights": (0.28, 0.01, 0.71, 1.0)},
    "saccharide": {"weights": (0.0, 1.0, 0.0, 1.0)},
    "terpene": {"weights": (0.2, 0.05, 0.75, 2.0)},
    "PKS-NRP_Hybrids": {"weights": (0.0, 0.22, 0.78, 1.0)},
    "other": {"weights": (0.01, 0.02, 0.97, 4.0)},
    "mix": {"weights": (0.2, 0.05, 0.75, 2.0)},
}


class RecordPairGenerator:
    """Generator to generate all-vs-all comparisons form a list of BGC records

    Attributes:
        label (str): Label for this bin
        source_records (list[BGCRecord]): List of BGC records to generate pairs from
    """

    def __init__(
        self,
        label: str,
        edge_param_id: int,
        weights: Optional[str] = None,
        record_type: Optional[RECORD_TYPE] = None,
    ):
        self.label = label
        self.edge_param_id = edge_param_id
        self.source_records: list[BGCRecord] = []
        self.record_ids: set[int] = set()
        if weights is None:
            weights = label
        self.weights = weights
        self.record_type = record_type

    def generate_pairs(
        self, legacy_sorting=False
    ) -> Generator[tuple[BGCRecord, BGCRecord], None, None]:
        """Returns a generator for all vs all record pairs in this bins

        This will always generate all pairs, and does not take into account any edges
        that already exist in the database

        Args:
            legacy_sorting (bool, optional): Whether to sort the BGC records by GBK file name.
            This is done in BiG-SCAPE 1.0 and can affect scoring depending on which of
            the BGC records is record A in a pair. TODO: may be removed in the future

        Yields:
            Generator[tuple[int, int]]: Generator for record pairs in this bin
        """
        for record_a, record_b in combinations(self.source_records, 2):
            if record_a.parent_gbk == record_b.parent_gbk:
                continue
            if legacy_sorting:
                sorted_a, sorted_b = sorted((record_a, record_b), key=sort_name_key)
                if sorted_a._db_id is None or sorted_b._db_id is None:
                    raise RuntimeError("generated pair is missing DB ids!")
                pair = (sorted_a, sorted_b)

            else:
                if record_a._db_id is None or record_b._db_id is None:
                    raise RuntimeError("generated pair is missing DB ids!")
                pair = (record_a, record_b)

            yield pair

    def generate_pair_ids(self, legacy_sorting=False):
        for pair in self.generate_pairs(legacy_sorting):
            yield (pair[0]._db_id, pair[1]._db_id)

    def num_pairs(self) -> int:
        """Return the number of pairs expected to be generated by the pairs Generator

        Returns:
            int: The number of pairs expected to be generated from the Generator
        """

        if len(self.source_records) < 2:
            return 0

        len_all_records = len(self.source_records)

        # (n*(n-1)) / 2
        num_all_pairs = int((len_all_records * (len_all_records - 1)) / 2)

        if self.record_type is not None and self.record_type != RECORD_TYPE.REGION:
            # if in protocluster mode this will overestimate the number of pairs as two
            # records from the same gbk will not be compared -> will not be a pair. Use
            # the database to find how many subrecords come from the same genbank, i.e.
            # how many pairs should be removed
            if not DB.metadata:
                raise RuntimeError("DB metadata is None!")
            record_table = DB.metadata.tables["bgc_record"]

            temp_record_id_table = create_temp_record_id_table(self.record_ids)

            # find a collection of gbks with more than one subrecord
            member_table = (
                select(func.count(record_table.c.gbk_id).label("rec_count"))
                .where(record_table.c.id.in_(select(temp_record_id_table.c.record_id)))
                .group_by(record_table.c.gbk_id)
                .having(func.count() > 1)
                .subquery()
            )
            # count how many times each occurs, e.g. gbks with 3 records occur 4 times
            subrecord_membership = select(
                member_table.c.rec_count, func.count(member_table.c.rec_count)
            ).group_by(member_table.c.rec_count)

            subrecord_counts = DB.execute(subrecord_membership).fetchall()

            # for each occurence of a number of subrecords, remove pairs
            if subrecord_counts is not None:
                for row in subrecord_counts:
                    nr_subrec, occurence = row
                    num_all_pairs -= int((nr_subrec * (nr_subrec - 1)) / 2 * occurence)

        return num_all_pairs

    def add_records(self, record_list: list[BGCRecord]):
        """Adds BGC records to this bin and creates a generator for the pairs

        Args:
            record_list (list[BGCRecord]): List of BGC records to add to this bin
        """
        self.source_records.extend(record_list)
        self.record_ids.update([record._db_id or -1 for record in record_list])

        # throw a ValueError if any region db id is None, as we expect all records to be
        # represented in the database
        if None in self.record_ids:
            raise ValueError("Record in bin has no db id!")

    def cull_singletons(self, cutoff: float, run_id: int):
        """Culls singletons for given cutoff, i.e. records that do not occur in any
        connected components

        Args:
            cutoff (float): distance cuttoff
            run_id (int): id of the current run

        Raises:
            RuntimeError: DB.metadata is None
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        cc_table = DB.metadata.tables["connected_component"]

        # get all record ids that occur in a connected component
        select_statement = (
            select(cc_table.c.record_id)
            .where(cc_table.c.cutoff == cutoff)
            .where(cc_table.c.bin_label == self.label)
            .where(cc_table.c.run_id == run_id)
        )

        connected_records = set(DB.execute(select_statement).scalars())

        self.record_ids = connected_records
        self.source_records = [
            record
            for record in self.source_records
            if record._db_id in connected_records
        ]

    def get_query_source_record_ids(self) -> list[int]:
        """Return a list of record ids of all QUERY source type records in this bin

        Returns:
            list[int]: record database ids
        """
        return [
            record._db_id
            for record in self.source_records
            if record._db_id
            and record.parent_gbk
            and record.parent_gbk.source_type == SOURCE_TYPE.QUERY
        ]

    def __repr__(self) -> str:
        return (
            f"Bin '{self.label}': {self.num_pairs()} pairs from "
            f"{len(self.source_records)} BGC records"
        )


class ConnectedComponentPairGenerator(RecordPairGenerator):
    """Generator that takes as input a conected component and generates
    the pairs from the edges in the component"""

    def __init__(self, connected_component, label: str):
        # getting the first one, assume consistent edge param id for all cc
        edge_param_id = connected_component[0][6]
        weights = bs_comparison.get_edge_weight(edge_param_id)

        super().__init__(label, edge_param_id, weights)
        self.connected_component = connected_component
        self.record_id_to_obj: dict[int, BGCRecord] = {}

    def add_records(self, record_list: list[BGCRecord]):
        """Adds BGC records to this bin and creates a generator for the pairs

        also creates a dictionary of record id to record objects
        """
        cc_record_ids = set()
        cc_record_list = []

        for edge in self.connected_component:
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

        for record in record_list:
            if record._db_id is None:
                raise ValueError("Region in bin has no db id!")
            if record._db_id not in cc_record_ids:
                continue

            self.record_id_to_obj[record._db_id] = record
            cc_record_list.append(record)

        return super().add_records(cc_record_list)

    def generate_pairs(
        self, legacy_sorting=False
    ) -> Generator[tuple[BGCRecord, BGCRecord], None, None]:
        """Returns an Generator for record pairs in this bin

        Args:
            legacy_sorting (bool, optional): Whether to sort the BGC records by GBK file name.
            This is done in BiG-SCAPE 1.0 and can affect scoring depending on which of
            the BGC records is record A in a pair.

        Yields:
            Generator[tuple[int, int]]: Generator for record pairs in this bin
        """

        for edge in self.connected_component:
            record_a_id, record_b_id, dist, jacc, adj, dss, edge_param_id = edge

            record_a = self.record_id_to_obj[record_a_id]
            record_b = self.record_id_to_obj[record_b_id]

            if record_a.parent_gbk == record_b.parent_gbk:
                continue

            if legacy_sorting:
                sorted_a, sorted_b = sorted((record_a, record_b), key=sort_name_key)
                if sorted_a._db_id is None or sorted_b._db_id is None:
                    raise RuntimeError("generated pair is missing DB ids!")
                pair = (sorted_a, sorted_b)
            else:
                if record_a._db_id is None or record_b._db_id is None:
                    raise RuntimeError("generated pair is missing DB ids!")
                pair = (record_a, record_b)

            yield pair


class MissingRecordPairGenerator(RecordPairGenerator):
    """Generator that wraps around another RecordPairGenerator to exclude any distances
    already in the database
    """

    def __init__(self, pair_generator: RecordPairGenerator):
        super().__init__(
            pair_generator.label, pair_generator.edge_param_id, pair_generator.weights
        )
        self.bin = pair_generator

    def num_pairs(self) -> int:
        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]

        temp_record_id_table = create_temp_record_id_table(self.bin.record_ids)

        # get all region._db_id in the bin where the record_a_id and record_b_id are in the
        # bin
        select_statement = (
            select(func.count(distance_table.c.record_a_id))
            .where(
                distance_table.c.record_a_id.in_(
                    select(temp_record_id_table.c.record_id)
                )
            )
            .where(
                distance_table.c.record_b_id.in_(
                    select(temp_record_id_table.c.record_id)
                )
            )
            .where(distance_table.c.edge_param_id == self.bin.edge_param_id)
        )

        # get count
        existing_distance_count = DB.execute(select_statement).scalar_one()

        # subtract from expected number of distances
        return self.bin.num_pairs() - existing_distance_count

    def generate_pairs(
        self, legacy_sorting=False
    ) -> Generator[tuple[BGCRecord, BGCRecord], None, None]:
        """Returns an Generator for record pairs in this bin

        Args:
            legacy_sorting (bool, optional): Whether to sort the BGC records by GBK file name.
            This is done in BiG-SCAPE 1.0 and can affect scoring depending on which of
            the BGC records is record A in a pair.

        Yields:
            Generator[tuple[int, int]]: Generator for record pairs in this bin
        """
        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]

        record_id_temp_table = create_temp_record_id_table(self.bin.record_ids)

        # get all region._db_id in the bin
        select_statement = (
            select(distance_table.c.record_a_id, distance_table.c.record_b_id)
            .where(
                distance_table.c.record_a_id.in_(
                    select(record_id_temp_table.c.record_id)
                )
            )
            .where(
                distance_table.c.record_b_id.in_(
                    select(record_id_temp_table.c.record_id)
                )
            )
            .where(distance_table.c.edge_param_id == self.bin.edge_param_id)
        )

        # generate a set of tuples of region id pairs
        existing_distances = set(DB.execute(select_statement).fetchall())

        for pair in self.bin.generate_pairs(legacy_sorting):
            # if the pair is not in the set of existing distances, yield it
            pair_ids = (pair[0]._db_id, pair[1]._db_id)
            if (
                pair_ids not in existing_distances
                and pair_ids[::-1] not in existing_distances
            ):
                yield pair

    def add_records(self, _: list[BGCRecord]):
        raise NotImplementedError("Cannot add records to a PartialRecordPairGenerator")


class QueryRecordPairGenerator(RecordPairGenerator):
    """Generator that generates pairs iteratively for query mode"""

    def __init__(
        self,
        label: str,
        edge_param_id: int,
        weights: str,
        record_type: Optional[RECORD_TYPE],
    ):
        super().__init__(label, edge_param_id, weights, record_type)
        self.reference_records: set[BGCRecord] = set()
        self.done_records: set[BGCRecord] = set()
        self.working_query_records: set[BGCRecord] = set()
        self.working_ref_records: set[BGCRecord] = set()
        self.record_id_to_obj: dict[int, BGCRecord] = {}

    def generate_pairs(
        self, legacy_sorting=False
    ) -> Generator[tuple[BGCRecord, BGCRecord], None, None]:
        """Returns a Generator for record pairs in this bin, pairs are generated in
        all vs all manner between working nodes and all reference nodes that are not
        done"""

        for record_a in self.working_query_records:
            for record_b in self.working_ref_records:
                if record_a == record_b:
                    continue

                if record_a.parent_gbk == record_b.parent_gbk:
                    continue

                if legacy_sorting:
                    sorted_a, sorted_b = sorted((record_a, record_b), key=sort_name_key)
                    if sorted_a._db_id is None or sorted_b._db_id is None:
                        raise RuntimeError("generated pair is missing DB ids!")
                    pair = (sorted_a, sorted_b)
                else:
                    if record_a._db_id is None or record_b._db_id is None:
                        raise RuntimeError("generated pair is missing DB ids!")
                    pair = (record_a, record_b)

                yield pair

    def num_pairs(self) -> int:
        """Returns the number of pairs expected to be generated by the pairs Generator,
        paris are generated in all vs all manner between working nodes and all reference
        nodes that are not done

        Returns:
            int: The number of pairs expected to be generated from the Generator
        """

        num_query_records = len(self.working_query_records)
        num_ref_records = len(self.working_ref_records)

        if len(self.source_records) < 2:
            return 0

        num_pairs = num_query_records * num_ref_records

        # delete pairs originating from the same parent gbk
        if self.record_type is not None and self.record_type != RECORD_TYPE.REGION:
            query_ids = [record._db_id for record in self.working_query_records]
            ref_ids = [record._db_id for record in self.working_ref_records]

            if DB.metadata is None:
                raise RuntimeError("DB.metadata is None")

            rec_table = DB.metadata.tables["bgc_record"]

            temp_record_id_table_query = create_temp_record_id_table(query_ids)
            temp_record_id_table_ref = create_temp_record_id_table(ref_ids)

            # contruct two tables that hold the gbk id and the number of subrecords
            # present in the set of query and ref records respectively
            query_gbk = (
                select(
                    rec_table.c.gbk_id,
                    func.count(rec_table.c.gbk_id).label("query_count"),
                )
                .where(
                    rec_table.c.id.in_(select(temp_record_id_table_query.c.record_id))
                )
                .group_by(rec_table.c.gbk_id)
                .subquery()
            )

            ref_gbk = (
                select(
                    rec_table.c.gbk_id,
                    func.count(rec_table.c.gbk_id).label("ref_count"),
                )
                .where(rec_table.c.id.in_(select(temp_record_id_table_ref.c.record_id)))
                .group_by(rec_table.c.gbk_id)
                .subquery()
            )

            # now we can join the two tables and obtain the number of links between
            # records from the same gbks by multiplying their counts
            same_gbk_query = select(
                func.sum(query_gbk.c.query_count * ref_gbk.c.ref_count)
            ).join(ref_gbk, query_gbk.c.gbk_id == ref_gbk.c.gbk_id)

            same_gbks = DB.execute(same_gbk_query).scalar_one()
            if same_gbks:
                num_pairs -= same_gbks

        return num_pairs

    def add_records(self, record_list: list[BGCRecord]) -> None:
        """Adds BGC records to this bin, additionaly splitting them into query and
        reference records

        also creates a dictionary of record id to record objects

        Args:
            record_list (list[BGCRecord]): List of BGC records to add to this bin
        """
        for record in record_list:
            if record._db_id is None:
                raise ValueError("Region in bin has no db id!")
            self.record_id_to_obj[record._db_id] = record

            if (
                record.parent_gbk is not None
                and record.parent_gbk.source_type == SOURCE_TYPE.QUERY
            ):
                self.working_query_records.add(record)
            else:
                self.reference_records.add(record)
                self.working_ref_records.add(record)

        super().add_records(record_list)

        return None

    def cycle_records(self, max_cutoff: float) -> None:
        """Cycle working records

        Resets the working sets of records for this bin, so that the working query
        records of the previous cycle now move to the done records, and for the pairs
        that have been just generated, those nodes in pairs that pass the max_cutoff distance
        threshold and arent in the done records are moved to the working query records,
        working ref records are also updated

        max_cutoff (float): highest cutoff used in the run
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]

        # TODO: in case these lists get too big for the sqlalchemy searches and that
        # crashes the run, add a temporary table to the database to store the working
        # query and ref records
        working_query_ids = [record._db_id for record in self.working_query_records]
        working_ref_ids = [record._db_id for record in self.working_ref_records]

        temp_record_id_table = create_temp_record_id_table(
            working_query_ids + working_ref_ids
        )

        select_statement = select(
            distance_table.c.record_a_id,
            distance_table.c.record_b_id,
            distance_table.c.distance,
        ).where(
            or_(
                and_(
                    distance_table.c.record_a_id.in_(
                        select(temp_record_id_table.c.record_id)
                    ),
                    distance_table.c.edge_param_id == self.edge_param_id,
                    distance_table.c.distance < max_cutoff,
                ),
                and_(
                    distance_table.c.record_a_id.in_(
                        select(temp_record_id_table.c.record_id)
                    ),
                    distance_table.c.edge_param_id == self.edge_param_id,
                    distance_table.c.distance < max_cutoff,
                ),
            )
        )

        # generate a set of tuples of region id pairs
        passing_distances = set(DB.execute(select_statement).fetchall())

        passing_records = set()

        for pair_ids in passing_distances:
            record_a = self.record_id_to_obj[pair_ids[0]]
            record_b = self.record_id_to_obj[pair_ids[1]]
            if record_a not in passing_records:
                passing_records.add(record_a)
            if record_b not in passing_records:
                passing_records.add(record_b)

        self.done_records.update(self.working_query_records)

        self.working_query_records = set(
            [record for record in passing_records if record not in self.done_records]
        )

        self.working_ref_records = set(
            [
                record
                for record in self.reference_records
                if record not in self.done_records
                and record not in self.working_query_records
            ]
        )

        return None


class QueryMissingRecordPairGenerator(RecordPairGenerator):
    """Generator that wraps around another RecordPairGenerator to exclude any distances
    already in the database"""

    def __init__(self, pair_generator: QueryRecordPairGenerator):
        super().__init__(
            pair_generator.label, pair_generator.edge_param_id, pair_generator.weights
        )
        self.bin = pair_generator

    def generate_pairs(
        self, legacy_sorting=False
    ) -> Generator[tuple[BGCRecord, BGCRecord], None, None]:
        """Returns an Generator for record pairs in this bin"""

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]

        working_query_ids = [record._db_id for record in self.bin.working_query_records]
        working_ref_ids = [record._db_id for record in self.bin.working_ref_records]

        # create a temporary table with the ids of the records in the bin
        temp_record_id_table = create_temp_record_id_table(
            working_query_ids + working_ref_ids
        )

        # get all region._db_id in the bin
        select_statement = select(
            distance_table.c.record_a_id, distance_table.c.record_b_id
        ).where(
            or_(
                and_(
                    distance_table.c.record_a_id.in_(select(temp_record_id_table)),
                    distance_table.c.edge_param_id == self.bin.edge_param_id,
                ),
                and_(
                    distance_table.c.record_a_id.in_(select(temp_record_id_table)),
                    distance_table.c.edge_param_id == self.bin.edge_param_id,
                ),
            )
        )

        # generate a set of tuples of region id pairs
        existing_distances = set(DB.execute(select_statement).fetchall())

        for pair in self.bin.generate_pairs(legacy_sorting):
            # if the pair is not in the set of existing distances, yield it
            pair_ids = (pair[0]._db_id, pair[1]._db_id)
            if (
                pair_ids not in existing_distances
                and pair_ids[::-1] not in existing_distances
            ):
                yield pair

    def num_pairs(self) -> int:
        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        distance_table = DB.metadata.tables["distance"]

        working_query_ids = [record._db_id for record in self.bin.working_query_records]
        working_ref_ids = [record._db_id for record in self.bin.working_ref_records]

        temp_record_id_table_query = create_temp_record_id_table(working_query_ids)
        temp_record_id_table_ref = create_temp_record_id_table(working_ref_ids)

        select_statement = select(func.count(distance_table.c.record_a_id)).where(
            or_(
                and_(
                    distance_table.c.record_a_id.in_(
                        select(temp_record_id_table_query.c.record_id)
                    ),
                    distance_table.c.record_b_id.in_(
                        select(temp_record_id_table_ref.c.record_id)
                    ),
                    distance_table.c.edge_param_id == self.bin.edge_param_id,
                ),
                and_(
                    distance_table.c.record_a_id.in_(
                        select(temp_record_id_table_ref.c.record_id)
                    ),
                    distance_table.c.record_b_id.in_(
                        select(temp_record_id_table_query.c.record_id)
                    ),
                    distance_table.c.edge_param_id == self.bin.edge_param_id,
                ),
            )
        )

        existing_distance_count = DB.execute(select_statement).scalar_one()

        return self.bin.num_pairs() - existing_distance_count

    def cycle_records(self, max_cutoff):
        self.bin.cycle_records(max_cutoff)


def generate_mix_bin(record_list: list[BGCRecord], run: dict) -> RecordPairGenerator:
    """Generate an all-vs-all bin of the supplied BGC records

    Args:
        bgc_list (list[BGCRecord]): BGC records to make into an all-vs-all bin

    Returns:
        BGCBin: The all-vs-all BGC bin
    """

    edge_param_id = bs_comparison.get_edge_param_id(run, "mix")

    record_type = run["record_type"]

    mix_bin = RecordPairGenerator(
        label="mix", edge_param_id=edge_param_id, record_type=record_type
    )

    mix_bin.add_records([record for record in record_list if record is not None])

    return mix_bin


def sort_name_key(record: BGCRecord) -> str:
    """Return the parent gbk file name without extension, or None if no parent gbk is
    assigned

    Args:
        record (BGCRecord): A BGCrecord

    Returns:
        str: the parent gbk file name without extension
    """
    if record.parent_gbk is None:
        return ""

    return record.parent_gbk.path.name[:-4]


def as_class_bin_generator(
    all_records: list[BGCRecord], run: dict
) -> Iterator[RecordPairGenerator]:
    """Generate bins for each antiSMASH class

    Args:
        gbks (list[GBK]): List of GBKs to generate bins for
        category_weights (str): weights to use for each class

    Yields:
        Iterator[RecordPairGenerator]: Generator that yields bins. Order is not guarenteed to be
        consistent
    """
    if run["legacy_weights"]:
        weight_type = "legacy_weights"
    else:
        weight_type = "mix"

    classify_mode = run["classify"]

    class_idx: dict[str, list[BGCRecord]] = {}
    category_weights: dict[str, str] = {}

    for record in all_records:
        # get region class for bin label and index
        if classify_mode == CLASSIFY_MODE.CLASS:
            record_class = record.product

        if classify_mode == CLASSIFY_MODE.CATEGORY:
            record_class = record.get_category()

        if run["hybrids_off"]:
            record_classes = record_class.split(".")
        else:
            record_classes = [record_class]

        for record_class in record_classes:
            if record_class not in class_idx:
                class_idx[record_class] = [record]

            if record_class in class_idx and record not in class_idx[record_class]:
                class_idx[record_class].append(record)

            if weight_type == "legacy_weights":
                # get region category for weights
                region_weight_cat = get_legacy_weights_from_category(
                    record, record_class, run
                )

                if record_class not in category_weights.keys():
                    category_weights[record_class] = region_weight_cat

            if weight_type == "mix":
                category_weights[record_class] = "mix"

    for class_name, records in class_idx.items():
        weight_category = category_weights[class_name]
        edge_param_id = bs_comparison.get_edge_param_id(run, weight_category)
        bin = RecordPairGenerator(
            class_name, edge_param_id, weight_category, run["record_type"]
        )
        bin.add_records(records)
        yield bin


def get_legacy_weights_from_category(
    record: BGCRecord, record_class: str, run: dict
) -> str:
    """Get the category of a BGC based on its antiSMASH product(s)
    and match it to the legacy weights classes

    Args:
        region (BGCRecord): region object

    Returns:
        str: class category to be used in weight selection
    """

    categories: list[str] = []

    if isinstance(record, ProtoCluster) or isinstance(record, ProtoCore):
        # T1PKS is the only case in which a antiSMASH category does not
        # correspond to a legacy_weights class
        if (
            record.category is not None
        ):  # for typing, we assume antismash 6 and up always have it
            if record.product == "T1PKS":
                categories.append(record.product)
            else:
                categories.append(record.category)

    if isinstance(record, Region):
        # get categories from region object
        for idx, cand_cluster in record.cand_clusters.items():
            if cand_cluster is not None:
                for idx, protocluster in cand_cluster.proto_clusters.items():
                    if protocluster is not None and protocluster.category is not None:
                        if protocluster.product == "T1PKS":
                            pc_category = protocluster.product
                        else:
                            pc_category = protocluster.category
                        # avoid duplicates, hybrids of the same kind use the same weight class
                        if run["hybrids_off"] and protocluster.product == record_class:
                            if pc_category not in categories:
                                categories.append(pc_category)
                        else:
                            if pc_category not in categories:
                                categories.append(pc_category)

    # for versions that dont have category information
    if len(categories) == 0:
        logging.warning(
            "No category found for %s",
            record,
            "This should not happen as long as antiSMASH is run with"
            "version 6 or up, consider whether there is something"
            "special about this region",
        )
        category = "other"

    # process into legacy_weights classes
    if len(categories) == 1:
        category = categories[0]

    if len(categories) > 1:
        if "NRPS" and ("PKS" or "T1PKS") in categories:
            category = "PKS-NRP_Hybrids"

        if "PKS" or ("PKS" and "T1PKS") in categories:
            category = "PKSother"  # PKS hybrids

        else:
            category = "other"  # other hybrids

    return category


def legacy_bin_generator(
    all_records: list[BGCRecord], run: dict
) -> Iterator[RecordPairGenerator]:  # pragma no cover
    """Generate bins for each class as they existed in the BiG-SCAPE 1.0 implementation

    Args:
        gbks (list[GBK]): List of GBKs to generate bins for

    Yields:
        Iterator[BGCBin]: Generator that yields bins. Order is not guarenteed to be
        consistent
    """
    # generate index
    class_idx: dict[str, list[BGCRecord]] = {}

    for record in all_records:
        if record is None:
            continue
        if record.product is None:
            continue

        product = record.product

        if run["hybrids_off"]:
            record_products = product.split(".")

        else:
            record_products = [product]

        for product in record_products:
            record_class = legacy_get_class(product)

            if record_class not in class_idx:
                class_idx[record_class] = [record]

            if record_class in class_idx and record not in class_idx[record_class]:
                class_idx[record_class].append(record)

    for class_name, records in class_idx.items():
        edge_param_id = bs_comparison.get_edge_param_id(run, class_name)
        bin = RecordPairGenerator(
            class_name, edge_param_id, class_name, run["record_type"]
        )
        bin.add_records(records)
        yield bin


# one of the few direct copy-and-pastes!
# TODO: test
def legacy_get_class(product):  # pragma no cover
    """Sort BGC by its type. Uses AntiSMASH annotations
    (see https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)

    Args:
        product (str): product type

    Returns:
        str: product class
    """

    # PKS_Type I
    if product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pks1_products"]:
        return "PKSI"
    # PKS Other Types
    elif product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pksother_products"]:
        return "PKSother"
    # NRPs
    elif product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["nrps_products"]:
        return "NRPS"
    # RiPPs
    elif product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["ripps_products"]:
        return "RiPP"
    # Saccharides
    elif product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["saccharide_products"]:
        return "saccharide"
    # Terpenes
    elif product == "terpene":
        return "terpene"
    # PKS/NRP hybrids
    elif len(product.split(".")) > 1:
        # print("  Possible hybrid: (" + cluster + "): " + product)
        # cf_fatty_acid category contains a trailing empty space

        subtypes = set(s.strip() for s in product.split("."))
        if (
            len(
                subtypes
                - (
                    BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pks1_products"]
                    | BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pksother_products"]
                    | BigscapeConfig.LEGACY_ANTISMASH_CLASSES["nrps_products"]
                )
            )
            == 0
        ):
            if (
                len(subtypes - BigscapeConfig.LEGACY_ANTISMASH_CLASSES["nrps_products"])
                == 0
            ):
                return "NRPS"
            elif (
                len(
                    subtypes
                    - (
                        BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pks1_products"]
                        | BigscapeConfig.LEGACY_ANTISMASH_CLASSES["pksother_products"]
                    )
                )
                == 0
            ):
                return "PKSother"  # pks hybrids
            else:
                return "PKS-NRP_Hybrids"
        elif (
            len(subtypes - BigscapeConfig.LEGACY_ANTISMASH_CLASSES["ripps_products"])
            == 0
        ):
            return "RiPP"
        elif (
            len(
                subtypes
                - BigscapeConfig.LEGACY_ANTISMASH_CLASSES["saccharide_products"]
            )
            == 0
        ):
            return "saccharide"
        else:
            return "other"  # other hybrid
    # Others
    elif product in BigscapeConfig.LEGACY_ANTISMASH_CLASSES["others_products"]:
        return "other"
    # ??
    elif product == "":
        # No product annotation. Perhaps not analyzed by antiSMASH
        return "other"
    else:
        logging.warning("unknown product %s", product)
        return "other"
