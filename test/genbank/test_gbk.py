"""Contains tests for the GBK class and functions"""
from pathlib import Path
from unittest import TestCase

from src.genbank.gbk import GBK


class TestGBK(TestCase):
    """Test class for base GBK parsing tests"""

    def test_parse_gbk(self):
        """Tests whether a GBK is instantiated correctly"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse_gbk(gbk_file_path)

        self.assertIsInstance(gbk, GBK)
