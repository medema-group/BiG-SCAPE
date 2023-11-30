"""module containing a test for loading files"""
# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.enums import SOURCE_TYPE
import big_scape.enums as bs_enums

# from this module
from big_scape.file_input.load_files import load_dataset_folder, load_gbk


class TestLoadGBK(TestCase):
    """Test class for loading of genbank files"""

    def test_gbk_path_invalid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a file, not a folder

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/valid_input.gbk"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
        }

        with self.assertRaises(NotADirectoryError):
            load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

    def test_gbk_path_valid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a folder containing a valid gbk file

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 4
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_gbk_path_empty(self):
        """Tests whether loading a given path returns none if the folder is empty"""
        # path pointing to a folder containing no gbk files

        run = {
            "input_dir": Path("test/test_data/empty_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        with self.assertRaises(FileNotFoundError):
            load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

    def test_load_gbk_not_a_file(self):
        """Tests whether the load_gbk function correctly returns none when input is not a file"""
        # path pointing to a folder, not a file
        gbk_file_path = Path("test/test_data/valid_gbk_folder/")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }
        with self.assertRaises(IsADirectoryError):
            load_gbk(gbk_file_path, SOURCE_TYPE.QUERY, run)

    def test_load_gbk_valid(self):
        """Tests whether the load_gbk function correctly loads a valid file

        nb: this does not test the contents of the GBK file, just the existence and reading of it
        """

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }
        gbk = load_gbk(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsNot(gbk, None)

    def test_include_exclude_str(self):
        """Tests whether include and exclude gbk strings filter correclty"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 1
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_include_all_override_exclude(self):
        """Tests whether include_gbk * correclty includes all files in dir"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["*"],
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 3
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)
