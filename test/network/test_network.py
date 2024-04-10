"""Contains tests for the networking module, which contains functions to retrieve
nodes and edges from the database
"""

# from python
import pathlib
from unittest import TestCase
from itertools import combinations

# from dependencies
from sqlalchemy import select, func

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.file_input as bs_files
import big_scape.data as bs_data
import big_scape.network.network as bs_network
import big_scape.comparison as bs_comparison
import big_scape.enums as bs_enums


def create_mock_gbk(i) -> bs_gbk.GBK:
    gbk = bs_gbk.GBK(
        pathlib.Path(f"test_path_{i}.gbk"), str(i), bs_enums.SOURCE_TYPE.QUERY
    )
    cds = bs_gbk.CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = bs_gbk.Region(gbk, i, 0, 100, False, "test")
    gbk.metadata = {"organism": "test", "taxonomy": "test", "description": "test"}
    return gbk


def gen_mock_edge_list(
    edge_gbks: list[bs_gbk.GBK],
    score: float = 0.0,
) -> list[
    tuple[int, int, float, float, float, float, int, bs_comparison.ComparableRegion]
]:
    edges = []
    for gbk_a, gbk_b in combinations(edge_gbks, 2):
        if gbk_a.region is None or gbk_b.region is None:
            continue
        if gbk_a.region._db_id is None or gbk_b.region._db_id is None:
            continue

        edges.append(
            (
                gbk_a.region._db_id,
                gbk_b.region._db_id,
                score,
                1.0 - score,
                1.0 - score,
                1.0 - score,
                1,
                bs_comparison.ComparableRegion(
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    0,
                    False,
                ),
            )
        )

    return edges


def create_mock_edge(a_id, b_id):
    return (
        a_id,
        b_id,
        0.0,
        1.0,
        1.0,
        1.0,
        1,
        bs_comparison.ComparableRegion(
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            False,
        ),
    )


class TestNetwork(TestCase):
    """Contains tests for the networking module"""

    def clear_db(self):
        """Clear the database"""
        bs_data.DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clear_db)

    def test_get_edge(self):
        """Test the get_edge function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(10):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # gather all database ids to include
        gbk_db_ids = set(gbk.region._db_id for gbk in gbks)

        # create a bunch of edges
        edges = gen_mock_edge_list(gbks)

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # providing an empty set just gets the first edge
        # which we expect to be (1, 2)
        expected_edge = (1, 2, 0.0, 1.0, 1.0, 1.0, 1)

        actual_edge = bs_network.get_edge(
            include_nodes=gbk_db_ids, edge_param_id=1, exclude_nodes=set()
        )

        self.assertEqual(expected_edge, actual_edge)

    def test_get_edge_with_exclusion(self):
        """Test the get_edge function with an exclusion"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(10):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # gather all database ids to include
        gbk_db_ids = set(gbk.region._db_id for gbk in gbks)

        # create a bunch of edges
        edges = gen_mock_edge_list(gbks)

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # we will ignore any edges that contain node with id 2
        # so the first one should be (1, 3)
        expected_edge = (1, 3, 0.0, 1.0, 1.0, 1.0, 1)

        actual_edge = bs_network.get_edge(
            include_nodes=gbk_db_ids, edge_param_id=1, exclude_nodes={2}
        )

        self.assertEqual(expected_edge, actual_edge)

    def test_get_edge_with_subselection(self):
        """Test the get_edge function with a subselection of included nodes"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(10):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # gather half of database ids to include
        gbk_db_ids = set(gbk.region._db_id for gbk in gbks[5:])

        # create a bunch of edges
        edges = gen_mock_edge_list(gbks)

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # we will not ignore any edges but only include half of gbks
        # so the first one should be (6, 7)
        expected_edge = (6, 7, 0.0, 1.0, 1.0, 1.0, 1)

        actual_edge = bs_network.get_edge(
            include_nodes=gbk_db_ids, edge_param_id=1, exclude_nodes={}
        )

        self.assertEqual(expected_edge, actual_edge)

    def test_get_edges(self):
        """Test the get_edges function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(10):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # create a bunch of edges
        edges = gen_mock_edge_list(gbks)

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # gen_mock_edge_list creates edges for all combinations of gbks
        # we'll include the first edge only. that should give us these edges:
        expected_edges = [
            (1, 2, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 3, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 4, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 5, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 6, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 7, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 8, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 9, 0.0, 1.0, 1.0, 1.0, 1),
            (1, 10, 0.0, 1.0, 1.0, 1.0, 1),
        ]

        actual_edges = bs_network.get_edges(set([1]))

        self.assertEqual(expected_edges, actual_edges)

    def test_generate_cc(self):
        """Test the generate_connected_components function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        # we will create 2 connected components

        gbks_a = []
        gbks_b = []
        # 0, 2, 4, 6, 8 are connected
        for i in range(5):
            gbk = create_mock_gbk(i * 2)
            gbks_a.append(gbk)
            gbk.save_all()

        # 1, 3, 5, 7, 9 are connected
        for i in range(5):
            gbk = create_mock_gbk(i * 2 + 1)
            gbks_b.append(gbk)
            gbk.save_all()

        # create edges for each connected component
        edges_a = gen_mock_edge_list(gbks_a)
        edges_b = gen_mock_edge_list(gbks_b)

        # save the edges
        for edge in edges_a:
            bs_comparison.save_edge_to_db(edge)
        for edge in edges_b:
            bs_comparison.save_edge_to_db(edge)

        # get a list of gbk ids to include

        # generate connected components
        bs_network.generate_connected_components(1, 1)

        # get the distinct connected components
        q = """
        SELECT COUNT(distinct id)
        FROM connected_component
        """

        expected_connected_components = 2
        actual_connected_components = bs_data.DB.execute_raw_query(q).fetchone()[0]

        self.assertEqual(expected_connected_components, actual_connected_components)

    def test_get_cc_ids(self):
        """Test the get_connected_component_ids function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        # we will create 2 connected components

        gbks_a = []
        a_ids = set()
        gbks_b = []
        b_ids = set()
        # 0, 2, 4, 6, 8 are connected
        for i in range(5):
            gbk = create_mock_gbk(i * 2)
            gbks_a.append(gbk)
            gbk.save_all()
            a_ids.add(gbk.region._db_id)

        # 1, 3, 5, 7, 9 are connected
        for i in range(5):
            gbk = create_mock_gbk(i * 2 + 1)
            gbks_b.append(gbk)
            gbk.save_all()
            b_ids.add(gbk.region._db_id)

        # create edges for each connected component
        edges_a = gen_mock_edge_list(gbks_a)
        edges_b = gen_mock_edge_list(gbks_b)

        # save the edges
        for edge in edges_a:
            bs_comparison.save_edge_to_db(edge)
        for edge in edges_b:
            bs_comparison.save_edge_to_db(edge)

        # get a list of gbk ids to include

        # generate connected components
        bs_network.generate_connected_components(1, 1)

        # get the cc ids
        cc_ids = bs_network.get_connected_component_ids(1, 1)

        # I insist on only having one assert per test
        # so we'll make things difficult for ourselves

        # these should both be 1
        num_a_in_cc = len(a_ids.intersection(cc_ids))
        num_b_in_cc = len(b_ids.intersection(cc_ids))

        expected_cc = 2

        # so this should be 2
        actual_cc = num_a_in_cc + num_b_in_cc

        self.assertEqual(expected_cc, actual_cc)

    def test_get_random_edge(self):
        """Test the get_random_edge function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        # we will create 2 connected components

        gbks_a = []
        ids = set()

        for i in range(5):
            gbk = create_mock_gbk(i)
            gbks_a.append(gbk)
            gbk.save_all()
            ids.add(gbk.region._db_id)

        # create edges for each connected component
        edges = gen_mock_edge_list(gbks_a)

        # these are full edges, reduce down to only the ids
        edges_small = [(edge[0], edge[1]) for edge in edges]

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        random_edge = bs_network.get_random_edge(1, 1)

        # we expect the edge to be in the list of edges
        self.assertIn(random_edge, edges_small)

    def test_create_temp_record_table(self):
        """Test the create_temp_record_table function"""

        # self.skipTest("Not implemented")

        bs_data.DB.create_in_mem()

        run = {"record_type": bs_enums.RECORD_TYPE.REGION}

        gbks = []

        for i in range(5):
            gbk = create_mock_gbk(i)
            gbk.save_all()
            gbks.append(gbk)

        include_records = []

        record = bs_files.get_all_bgc_records(run, gbks)
        include_records.extend(record)

        table = bs_network.create_temp_record_table(include_records)

        temp_table_name = table.name
        temp_table = bs_data.DB.metadata.tables[temp_table_name]

        rows = bs_data.DB.execute_raw_query(
            "select * from " + temp_table_name
        ).fetchall()
        row_count_query = select(func.count(temp_table.c.record_id))

        expected_row_count = 5

        actual_row_count = bs_data.DB.execute(row_count_query).fetchone()[0]

        self.assertEqual(expected_row_count, actual_row_count)
        self.assertEqual(expected_row_count, len(rows))

    def test_get_cc_edges(self):
        """Test the get_cc_edges function"""

        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks_a = []
        for i in range(3):
            gbk = create_mock_gbk(i)
            gbks_a.append(gbk)
            gbk.save_all()

        gbks_b = []
        for i in range(3):
            gbk = create_mock_gbk(i)
            gbks_b.append(gbk)
            gbk.save_all()

        # create a bunch of edges, generating a cc with edges of distance 0
        # and another with distance 0.6
        # gen_mock_edge_list creates edges for all combinations of gbks
        edges_a = gen_mock_edge_list(gbks_a)
        edges_b = gen_mock_edge_list(gbks_b, 0.6)

        edges = edges_a + edges_b

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # generate connected components
        bs_network.generate_connected_components(0.5, 1)

        cc_ids = bs_network.get_connected_component_ids(0.5, 1)

        cc_edges_a = bs_network.get_cc_edges(cc_ids[0], 1, 1)

        expected_data = {"cc_count": 1, "edge_count": len(edges_a)}
        seen_data = {"cc_count": len(cc_ids), "edge_count": len(cc_edges_a)}

        self.assertEqual(expected_data, seen_data)

    def test_get_nodes_from_cc(self):
        """Test the get_nodes_from_cc function"""

        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks_a = []
        for i in range(3):
            gbk = create_mock_gbk(i)
            gbks_a.append(gbk)
            gbk.save_all()

        gbks_b = []
        for i in range(3):
            gbk = create_mock_gbk(i)
            gbks_b.append(gbk)
            gbk.save_all()

        # gen_mock_edge_list creates edges for all combinations of gbks
        edges_a = gen_mock_edge_list(gbks_a)

        # save the edges
        for edge in edges_a:
            bs_comparison.save_edge_to_db(edge)

        # generate connected components
        cc = []
        for edge in edges_a:
            (
                record_a_id,
                record_b_id,
                dist,
                jacc,
                adj,
                dss,
                edge_param_id,
                comparable_region,
            ) = edge

            shortened = (record_a_id, record_b_id, dist, jacc, adj, dss, edge_param_id)
            cc.append(shortened)

        cc_nodes = bs_network.get_nodes_from_cc(cc, gbks_a + gbks_b)

        self.assertEqual(len(cc_nodes), len(gbks_a))
