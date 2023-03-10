"""Contains tests for the GBK class and functions"""
from pathlib import Path
from unittest import TestCase

from src.genbank.gbk import GBK
from src.genbank.region import Region


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

        self.assertRaises(ValueError, GBK.parse, gbk_file_path)

    def test_populate_regions(self):
        """Tests whether parsing a GBK correctly populates the underlying regions"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        self.assertIsInstance(gbk.region, Region)
