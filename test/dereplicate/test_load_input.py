"""module containing unit tests for loading input files"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.SeqRecord import SeqRecord

# from other modules
import big_scape.enums as bs_enums

# from this module
from big_scape.dereplicating.data_loading import load_input_folder, parse_gbk_files


class TestLoadInput(TestCase):
    """Test class for loading files and records in files"""

    def test_load_input_folder(self):
        """Tests whether the load_input_folder function correctly loads and filters files"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
        }

        input_gbks = load_input_folder(run)

        self.assertEqual(len(input_gbks), 2)

    def test_parse_gbk_files(self):
        """Tests whether the parse_gbk_files function correctly parses GBK files"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
        }

        input_gbk_paths = load_input_folder(run)

        gbk_data = parse_gbk_files(input_gbk_paths)

        gbk_data_list = list(gbk_data)

        self.assertIsInstance(gbk_data_list[0], tuple)
        self.assertIsInstance(gbk_data_list[0][0], Path)
        self.assertIsInstance(gbk_data_list[0][1], str)
        self.assertIsInstance(gbk_data_list[0][2], SeqRecord)

