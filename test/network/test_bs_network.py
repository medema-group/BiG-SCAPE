"""Contains tests for bigscape networking classes and methods"""

# from python
from pathlib import Path
from unittest import TestCase
from itertools import combinations
from src.data.sqlite import DB

# from other modules
from src.genbank import GBK, Region, BGCRecord
from src.comparison import BGCPair
from src.network import BSNetwork
from src.enums import SOURCE_TYPE


class TestBSNetwork(TestCase):
    """Contains tests for the BSNetwork class"""

    def test_add_node(self):
        """Tests whether a node can be created from a BGCRecord object"""
        gbk = GBK("", "test")
        region = Region(gbk, 0, 0, 100, False, "")

        network = BSNetwork()

        network.add_node(region)

        self.assertEqual(network.graph.number_of_nodes(), 1)

    def test_get_nodes_no_source(self):
        """Tests whether the correct nodes are returned when get_nodes is called"""

        parent_gbk_query = GBK(Path("test"), source_type=SOURCE_TYPE.QUERY)
        parent_gbk_ref = GBK(Path("test"), source_type=SOURCE_TYPE.REFERENCE)
        parent_gbk_mibig = GBK(Path("test"), source_type=SOURCE_TYPE.MIBIG)
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_mibig, 0, 0, 10, False, "")

        network = BSNetwork()
        network.add_node(bgc_a)
        network.add_node(bgc_b)
        network.add_node(bgc_c)
        network.add_node(bgc_d)

        expected_nodes = [bgc_a, bgc_b, bgc_c, bgc_d]
        nodes = network.get_nodes()

        self.assertEqual(expected_nodes, nodes)

    def test_get_nodes_with_source(self):
        """Tests whether the correct nodes are returned"""

        parent_gbk_query = GBK(Path("test"), source_type=SOURCE_TYPE.QUERY)
        parent_gbk_ref = GBK(Path("test"), source_type=SOURCE_TYPE.REFERENCE)
        parent_gbk_mibig = GBK(Path("test"), source_type=SOURCE_TYPE.MIBIG)
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_mibig, 0, 0, 10, False, "")

        network = BSNetwork()
        network.add_node(bgc_a)
        network.add_node(bgc_b)
        network.add_node(bgc_c)
        network.add_node(bgc_d)

        expected_nodes = [bgc_a, bgc_b, bgc_d]
        nodes = network.get_nodes(node_types=[SOURCE_TYPE.QUERY, SOURCE_TYPE.MIBIG])

        self.assertEqual(expected_nodes, nodes)

    def test_get_singleton_nodes_with_source(self):
        """Tests whether the correct singletons are returned"""

        parent_gbk_query = GBK(Path("test"), source_type=SOURCE_TYPE.QUERY)
        parent_gbk_ref = GBK(Path("test"), source_type=SOURCE_TYPE.REFERENCE)
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")

        network = BSNetwork()
        network.add_node(bgc_a)
        network.add_node(bgc_b)
        network.add_node(bgc_c)
        network.add_node(bgc_d)

        network.add_edge_pair(BGCPair(bgc_a, bgc_b))  # query to query
        network.add_edge_pair(BGCPair(bgc_a, bgc_c))  # query to ref

        expected_ref_singletons = [bgc_d]
        ref_singletons = network.get_singletons(
            node_types=[SOURCE_TYPE.REFERENCE, SOURCE_TYPE.MIBIG]
        )

        self.assertEqual(expected_ref_singletons, ref_singletons)

    def test_cull_singletons_no_source(self):
        """Tests whether all singletons are culled"""

        parent_gbk_query = GBK(Path("test"), source_type=SOURCE_TYPE.QUERY)
        parent_gbk_ref = GBK(Path("test"), source_type=SOURCE_TYPE.REFERENCE)
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")

        network = BSNetwork()
        network.add_node(bgc_a)
        network.add_node(bgc_b)
        network.add_node(bgc_c)
        network.add_node(bgc_d)

        network.cull_singletons()
        post_cull_nodes = network.get_nodes()
        expected_nodes = []

        self.assertEqual(post_cull_nodes, expected_nodes)

    def test_cull_singletons_with_source(self):
        """Tests whether the correct singletons are culled"""

        parent_gbk_query = GBK(Path("test"), source_type=SOURCE_TYPE.QUERY)
        parent_gbk_ref = GBK(Path("test"), source_type=SOURCE_TYPE.REFERENCE)
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")

        network = BSNetwork()
        network.add_node(bgc_a)
        network.add_node(bgc_b)
        network.add_node(bgc_c)
        network.add_node(bgc_d)

        network.cull_singletons(node_types=[SOURCE_TYPE.QUERY])
        post_cull_nodes = network.get_nodes()
        expected_nodes = [bgc_c, bgc_d]

        self.assertEqual(post_cull_nodes, expected_nodes)

    def test_add_edge(self):
        """Tests whether an edge can be correctly added between two nodes"""
        gbk_a = GBK("", "test")
        region_a = Region(gbk_a, 0, 0, 100, False, "")

        gbk_b = GBK("", "test")
        region_b = Region(gbk_b, 0, 0, 100, False, "")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge_pair(pair)

        self.assertEqual(network.graph.number_of_edges(), 1)

    def test_add_edge_nonexistent_node(self):
        """Tests whether add_edge correctly raises a ValueError if one node is not yet
        present in the graph
        """
        gbk_a = GBK(Path("test/test_data/tmp/a.gbk"), "test")
        region_a = Region(gbk_a, 0, 0, 100, False, "")

        gbk_b = GBK(Path("test/test_data/tmp/b.gbk"), "test")
        region_b = Region(gbk_b, 0, 0, 100, False, "")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        # only add a
        network.add_node(region_a)

        self.assertRaises(KeyError, network.add_edge_pair, pair)

    def test_write_graphml(self):
        """Tests whether the graph can be written to an output grpahm file"""

        gbk_a = GBK(Path("test1.gbk"), "test")
        region_a = Region(gbk_a, 0, 0, 100, False, "")

        gbk_b = GBK(Path("test2.gbk"), "test")
        region_b = Region(gbk_b, 0, 0, 100, False, "")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge_pair(pair)

        test_out_path = Path("test/test_data/tmp/test.graphml")

        network.write_graphml(test_out_path)

        self.assertTrue(test_out_path.exists())

    def test_write_tsv(self):
        """Tests whether the graph can be written to an output grpahm file"""

        gbk_a = GBK(Path("test1.gbk"), "test")
        region_a = Region(gbk_a, 0, 0, 100, False, "")

        gbk_b = GBK(Path("test2.gbk"), "test")
        region_b = Region(gbk_b, 0, 0, 100, False, "")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge_pair(pair)

        test_out_path = Path("test/test_data/tmp/test.tsv")

        network.write_edgelist_tsv(test_out_path)

        self.assertTrue(test_out_path.exists())

    def test_gen_cutoff_subgraphs(self):
        """Tests whether the correct connected component subgraphs are created when
        a graph is given a cutoff"""

        network = BSNetwork()

        gbk = GBK(Path("test1.gbk"), "test")

        # 6 regions
        regions = []
        for region_number in range(6):
            region = Region(gbk, region_number, 0, 100, False, "")
            regions.append(region)
            network.add_node(region)

        edges = list(combinations(regions, 2))

        test_cutoff = 0.3

        # subgraphs will look like this:
        #    0     3
        #   / \   / \
        #  1---2 4---5
        subgraph_a = [regions[i] for i in [0, 1, 2]]
        subgraph_b = [regions[i] for i in [3, 4, 5]]
        for idx, edge in enumerate(edges):
            bgc_a, bgc_b = edge
            pair = BGCPair(bgc_a, bgc_b)

            if bgc_a in subgraph_a and bgc_b in subgraph_a:
                distance = 0.2
            elif bgc_a in subgraph_b and bgc_b in subgraph_b:
                distance = 0.2
            else:
                distance = 0.5

            network.add_edge_pair(pair, dist=distance)

        subgraphs = network.generate_cutoff_subgraphs("dist", test_cutoff)

        expected_subgraph_sizes = [3, 3]
        actual_subgraph_sizes = [len(subgraph) for subgraph in subgraphs]

        self.assertEqual(expected_subgraph_sizes, actual_subgraph_sizes)

    def test_load_from_db(self):
        """Tests whether a network can be correctly loaded from the database"""
        # open in-memory db
        DB.create_in_mem()

        network = BSNetwork()

        gbk = GBK(Path("test1.gbk"), "test")
        gbk.save()

        # 6 regions
        regions = []
        for region_number in range(6):
            region = Region(gbk, region_number, 0, 100, False, "")
            region.save()
            regions.append(region)
            network.add_node(region)

        edges_a = list(combinations(regions[:3], 2))
        edges_b = list(combinations(regions[3:], 2))

        # network will look like this:
        #   0     3
        #  / \   / \
        # 1---2 4---5
        #   A     B

        for edge in edges_a:
            pair = BGCPair(edge[0], edge[1])
            network.add_edge_pair(pair, dist=0.0, jc=0.0, ai=0.0, dss=0.0)

        for edge in edges_b:
            pair = BGCPair(edge[0], edge[1])
            network.add_edge_pair(pair, dist=0.0, jc=0.0, ai=0.0, dss=0.0)

        # save network to do
        network.export_distances_to_db()

        # load the network from db
        loaded_network = BSNetwork.load_from_db(set(regions))

        # check that the loaded network is the same as the original
        self.assertEqual(network, loaded_network)
