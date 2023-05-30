"""Contains tests for bigscape networking classes and methods"""

# from python
from pathlib import Path
from unittest import TestCase

# from other modules
from src.genbank import GBK, Region
from src.comparison import BGCPair
from src.network import BSNetwork


class TestBSNetwork(TestCase):
    """Contains tests for the BSNetwork class"""

    def test_add_node(self):
        """Tests whether a node can be created from a BGCRecord object"""
        region = Region(1)
        region.nt_start = 0
        region.nt_stop = 100
        region.parent_gbk = GBK("", "test")

        network = BSNetwork()

        network.add_node(region)

        self.assertEqual(network.graph.number_of_nodes(), 1)

    def test_add_edge(self):
        """Tests whether an edge can be correctly added between two nodes"""
        region_a = Region(1)
        region_a.nt_start = 0
        region_a.nt_stop = 100
        region_a.parent_gbk = GBK("", "test")

        region_b = Region(1)
        region_b.nt_start = 0
        region_b.nt_stop = 100
        region_b.parent_gbk = GBK("", "test")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge(pair)

        self.assertEqual(network.graph.number_of_edges(), 1)

    def test_add_edge_nonexistent_node(self):
        """Tests whether add_edge correctly raises a ValueError if one node is not yet
        present in the graph
        """
        region_a = Region(1)
        region_a.nt_start = 0
        region_a.nt_stop = 100
        region_a.parent_gbk = GBK(Path("test/test_data/tmp/a.gbk"), "test")

        region_b = Region(1)
        region_b.nt_start = 0
        region_b.nt_stop = 100
        region_b.parent_gbk = GBK(Path("test/test_data/tmp/b.gbk"), "test")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        # only add a
        network.add_node(region_a)

        self.assertRaises(KeyError, network.add_edge, pair)

    def test_write_graphml(self):
        """Tests whether the graph can be written to an output grpahm file"""

        region_a = Region(1)
        region_a.nt_start = 0
        region_a.nt_stop = 100
        region_a.product = ""
        region_a.contig_edge = False
        region_a.parent_gbk = GBK(Path("test1.gbk"), "test")

        region_b = Region(1)
        region_b.nt_start = 0
        region_b.nt_stop = 100
        region_b.product = ""
        region_b.contig_edge = False
        region_b.parent_gbk = GBK(Path("test2.gbk"), "test")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge(pair)

        test_out_path = Path("test/test_data/tmp/test.graphml")

        network.write_graphml(test_out_path)

        self.assertTrue(test_out_path.exists())

    def test_write_tsv(self):
        """Tests whether the graph can be written to an output grpahm file"""

        region_a = Region(1)
        region_a.nt_start = 0
        region_a.nt_stop = 100
        region_a.product = ""
        region_a.contig_edge = False
        region_a.parent_gbk = GBK(Path("test1.gbk"), "test")

        region_b = Region(1)
        region_b.nt_start = 0
        region_b.nt_stop = 100
        region_b.product = ""
        region_b.contig_edge = False
        region_b.parent_gbk = GBK(Path("test2.gbk"), "test")

        pair = BGCPair(region_a, region_b)

        network = BSNetwork()

        network.add_node(region_a)
        network.add_node(region_b)

        network.add_edge(pair)

        test_out_path = Path("test/test_data/tmp/test.tsv")

        network.write_edgelist_tsv(test_out_path)

        self.assertTrue(test_out_path.exists())
