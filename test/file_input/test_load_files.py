"""module containing a test for loading files"""
# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from src.genbank import SOURCE_TYPE

# from this module
from src.file_input.load_files import load_dataset_folder, load_gbk


class TestLoadGBK(TestCase):
    """Test class for loading of genbank files"""

    def test_gbk_path_invalid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a file, not a folder
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        with self.assertRaises(NotADirectoryError):
            load_dataset_folder(gbk_file_path, SOURCE_TYPE.QUERY)

    def test_gbk_path_valid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a folder containing a valid gbk file
        gbk_file_path = Path("test/test_data/valid_gbk_folder/")

        load_result = load_dataset_folder(gbk_file_path, SOURCE_TYPE.QUERY)

        expected_count = 3
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_gbk_path_empty(self):
        """Tests whether loading a given path returns none if the folder is empty"""
        # path pointing to a folder containing no gbk files
        gbk_file_path = Path("test/test_data/empty_gbk_folder/")

        with self.assertRaises(FileNotFoundError):
            load_dataset_folder(gbk_file_path, SOURCE_TYPE.QUERY)

    def test_load_gbk_not_a_file(self):
        """Tests whether the load_gbk function correctly returns none when input is not a file"""
        # path pointing to a folder, not a file
        gbk_file_path = Path("test/test_data/valid_gbk_folder/")

        with self.assertRaises(IsADirectoryError):
            load_gbk(gbk_file_path, SOURCE_TYPE.QUERY)

    def test_load_gbk_valid(self):
        """Tests whether the load_gbk function correctly loads a valid file

        nb: this does not test the contents of the GBK file, just the existence and reading of it
        """

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk = load_gbk(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsNot(gbk, None)

    def test_include_exclude_str(self):
        """Tests whether include and exclude gbk strings filter correclty"""

        gbk_folder_path = Path("test/test_data/alt_valid_gbk_input/")

        load_result = load_dataset_folder(gbk_folder_path, SOURCE_TYPE.QUERY)

        expected_count = 1
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_include_all_override_exclude(self):
        """Tests whether include_gbk * correclty includes all files in dir"""

        gbk_folder_path = Path("test/test_data/alt_valid_gbk_input/")

        load_result = load_dataset_folder(
            gbk_folder_path, SOURCE_TYPE.QUERY, include_gbk=["*"]
        )

        expected_count = 3
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)
