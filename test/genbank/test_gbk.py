"""Contains tests for the GBK class and functions"""
from pathlib import Path
from unittest import TestCase

from Bio.Seq import Seq

from src.genbank.gbk import GBK
from src.genbank.region import Region
from src.genbank.proto_core import Protocore
from src.errors.genbank import InvalidGBKError


class TestGBK(TestCase):
    """Test class for base GBK parsing tests"""

    def test_parse_gbk(self):
        """Tests whether a GBK is instantiated correctly"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        self.assertIsInstance(gbk, GBK)

    def test_parse_gbk_multiple_regions(self):
        """Tests whether a GBK file has more than one region"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_multiple_regions_folder/valid_input_multiple_regions.gbk"
        )

        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path)

    def test_populate_regions(self):
        """Tests whether parsing a GBK correctly populates the underlying region"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        self.assertIsInstance(gbk.region, Region)

    def test_populate_hierarchical_objects(self):
        """Tests whether parsing a GBK correclty generates parent-child feature relations
        via checking for presence of the lowest level child - proto_core"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        proto_core = gbk.region.cand_clusters[1].proto_clusters[1].protocore[1]

        self.assertIsInstance(proto_core, Protocore)

    def test_parse_gbk_has_dna_seq(self):
        """Tests whether parsing a GBK correclty has DNA sequence"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        dna_sequence = gbk.nt_seq

        self.assertIsInstance(dna_sequence, Seq)
