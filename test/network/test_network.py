"""Contains tests for the networking module, which contains functions to retrieve
nodes and edges from the database
"""

# from python
import pathlib
from unittest import TestCase
from itertools import combinations

# from other modules
import big_scape.genbank as bs_gbk
import big_scape.data as bs_data
import big_scape.network.network as bs_network
import big_scape.comparison as bs_comparison
import big_scape.enums as bs_enums


def create_mock_gbk(i) -> bs_gbk.GBK:
    gbk = bs_gbk.GBK(pathlib.Path(f"test_path_{i}.gbk"), bs_enums.SOURCE_TYPE.QUERY)
    cds = bs_gbk.CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = bs_gbk.Region(gbk, 1, 0, 100, False, "test")
    gbk.metadata = {"organism": "test", "taxonomy": "test", "description": "test"}
    return gbk


def gen_mock_edge_list(
    edge_gbks: list[bs_gbk.GBK],
) -> list[
    tuple[
        int,
        int,
        float,
        float,
        float,
        float,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        bool,
    ]
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
                0.0,
                1.0,
                1.0,
                1.0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                False,
            )
        )

    return edges


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

        actual_edge = bs_network.get_edge(include_nodes=gbk_db_ids, exclude_nodes=set())

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

        actual_edge = bs_network.get_edge(include_nodes=gbk_db_ids, exclude_nodes={2})

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

        actual_edge = bs_network.get_edge(include_nodes=gbk_db_ids, exclude_nodes={})

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

    def test_get_connected_edges(self):
        """Test the get_connected_edges function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(3):
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

        # take the first edge (1,2) as starting point
        seed_edge = bs_network.get_edge(include_nodes=gbk_db_ids, exclude_nodes=set())

        # nodes connected to starting edge to find new edges for: 1, 2
        connected_nodes = set(list(seed_edge[:2]))
        connected_component = set((seed_edge,))

        # find all edges connected to nodes 1 and 2, excluding starting edge (1,2)
        connected_edges = bs_network.get_connected_edges(
            include_nodes=gbk_db_ids,
            connected_nodes=connected_nodes,
            connected_component=connected_component,
        )

        # all edges connected to only nodes 1 and 2 are therefore (1,3) and (2,3)
        expected_edges = [(1, 3, 0.0, 1.0, 1.0, 1.0, 1), (2, 3, 0.0, 1.0, 1.0, 1.0, 1)]

        self.assertEqual(expected_edges, connected_edges)

    def test_get_connected_edges_with_subselection(self):
        """Test the get_connected_edges function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(4):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # gather database ids to include, exclude the last gbk entry
        gbk_db_ids = set(gbk.region._db_id for gbk in gbks[:-1])

        # create a bunch of edges
        edges = gen_mock_edge_list(gbks)

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # take the first edge (1,2) as starting point
        seed_edge = bs_network.get_edge(include_nodes=gbk_db_ids, exclude_nodes=set())

        # nodes connected to starting edge to find new edges for: 1, 2
        connected_nodes = set(list(seed_edge[:2]))
        connected_component = set((seed_edge,))

        # find all edges connected to nodes 1 and 2, excluding starting edge (1,2)
        # since gbk 4 is not in include_nodes, its connected edges are not returned
        connected_edges = bs_network.get_connected_edges(
            include_nodes=gbk_db_ids,
            connected_nodes=connected_nodes,
            connected_component=connected_component,
        )

        # all edges connected to only nodes 1 and 2 are therefore (1,3) and (2,3)
        expected_edges = [(1, 3, 0.0, 1.0, 1.0, 1.0, 1), (2, 3, 0.0, 1.0, 1.0, 1.0, 1)]

        self.assertEqual(expected_edges, connected_edges)

    def test_get_query_connected_component(self):
        """Test get_query_connected_component function"""
        bs_data.DB.create_in_mem()

        # create a bunch of gbk files
        gbks = []
        for i in range(6):
            gbk = create_mock_gbk(i)
            gbks.append(gbk)
            gbk.save_all()

        # gather all records to include
        gbk_regions = [gbk.region for gbk in gbks]

        # create two connected components
        edges = gen_mock_edge_list(gbks[:3])
        edges = gen_mock_edge_list(gbks[3:])

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # find cc for query record 5
        query_cc = bs_network.get_query_connected_component(
            include_records=gbk_regions, query_node_id=5
        )

        # only cc containing query should contain records 4, 5 and 6
        expected_cc = [
            (4, 5, 0.0, 1.0, 1.0, 1.0, 1),
            (4, 6, 0.0, 1.0, 1.0, 1.0, 1),
            (5, 6, 0.0, 1.0, 1.0, 1.0, 1),
        ]

        self.assertEqual(expected_cc, query_cc)
