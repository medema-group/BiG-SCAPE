"""Contains tests for parameter validation functions"""

# from python
from pathlib import Path
from unittest import TestCase
from random import randint
import time

# from other modules
from src.errors import InvalidArgumentError, ArgumentParseError
from src.parameters.input import (
    validate_input_dir,
    validate_input_mode,
    validate_cds_overlap_cutoff,
)
from src.parameters.hmmer import validate_includelist, validate_hsp_overlap_cutoff
from src.parameters.binning import validate_mix, validate_query_bgc
from src.parameters.comparison import validate_alignment_mode
from src.parameters.networking import validate_gcf_cutoffs
from src.parameters.output import (
    validate_output_dir,
    validate_db_path,
    validate_log_path,
)


class TestInputValidation(TestCase):
    """Contains tests on the input validation functions"""

    def test_validate_input_dir_none(self):
        """Tests whether validate_input_dir raises an exception when None is passed as
        path
        """
        self.assertRaises(InvalidArgumentError, validate_input_dir, None)

    def test_validate_input_dir_not_exists(self):
        """Tests whether validate_input_dir raises an exception when a path is passed
        that does not exist
        """
        missing_folder_path = Path("test/test_data/non_existent_path")

        self.assertRaises(InvalidArgumentError, validate_input_dir, missing_folder_path)

    def test_validate_input_dir_not_a_dir(self):
        """Tests whether validate_input_dir raises an exception when the location
        pointed to by path is not a directory
        """
        not_a_folder_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region.gbk"
        )
        self.assertRaises(InvalidArgumentError, validate_input_dir, not_a_folder_path)

    def test_validate_input_mode_wrong_mode(self):
        """Tests whether validate_input_mode raises an exception when the mode passed is
        not a valid input mode
        """
        wrong_mode = "pancakes"

        self.assertRaises(InvalidArgumentError, validate_input_mode, wrong_mode)

    def test_validate_cds_overlap_cutoff_low(self):
        """Tests whether validate_overlap_cutoff raises an exception if the cutoff is
        too low
        """
        self.assertRaises(InvalidArgumentError, validate_cds_overlap_cutoff, -1.0)

    def test_validate_cds_overlap_cutoff_high(self):
        """Tests whether validate_overlap_cutoff raises an exception if the cutoff is
        too low
        """
        self.assertRaises(InvalidArgumentError, validate_cds_overlap_cutoff, 2.0)

    def test_validate_query_bgc_path_not_exist(self):
        """Tests whether validate_query_bgc raises an exception if the given path does
        not exist"""
        missing_file_path = Path("test/test_data/non_existent_path")

        self.assertRaises(InvalidArgumentError, validate_query_bgc, missing_file_path)

    def test_validate_query_bgc_path_not_a_file(self):
        """Tests whether validate_query_bgc raises and exception if the given path is
        not a file"""
        not_a_file_path = Path("test/test_data/")

        self.assertRaises(InvalidArgumentError, validate_query_bgc, not_a_file_path)


class TestHMMerValidation(TestCase):
    """Contains tests on the HMMer validation functions"""

    def test_validate_includelist_valid_accessions(self):
        """Tests whether validate_domain_includelist correctly returns a list of strings
        when given a valid domain includelist file
        """
        valid_file_path = Path(
            "test/test_data/domain_includelist/valid_domain_includelist.txt"
        )

        expected_domains = [
            "PF00014",
            "PF00020.18",
        ]

        actual_domains = validate_includelist(valid_file_path)

        self.assertEqual(expected_domains, actual_domains)

    def test_validate_includelist_invalid_accession(self):
        """Tests whether validate_domain_includelist raises and exception if an invalid
        pfam domain accession is present in a given domain includelist file
        """
        invalid_file_path = Path(
            "test/test_data/domain_includelist/invalid_domain_includelist.txt"
        )

        self.assertRaises(ArgumentParseError, validate_includelist, invalid_file_path)

    def test_validate_hsp_overlap_cutoff_low(self):
        """Tests whether validate_overlap_cutoff raises an exception if the cutoff is
        too low
        """
        self.assertRaises(InvalidArgumentError, validate_hsp_overlap_cutoff, -1.0)

    def test_validate_hsp_overlap_cutoff_high(self):
        """Tests whether validate_overlap_cutoff raises an exception if the cutoff is
        too low
        """
        self.assertRaises(InvalidArgumentError, validate_hsp_overlap_cutoff, 2.0)


class TestBinningValidation(TestCase):
    """Contains tests on the binning validation functions"""

    def test_validate_mix_none(self):
        """Tests whether validate_mix raises an exception when mix is somehow set to
        none. This should never happen since the default value is set both in the class
        as well as in the argument parser
        """
        self.assertRaises(InvalidArgumentError, validate_mix, None)


class TestComparisonValidation(TestCase):
    """Contains tests on the comparison validation functions"""

    def test_validate_alignment_mode_wrong_mode(self):
        """Tests whether validate_alignment_mode raises an exception when the mode
        passed is not a valid comparison mode
        """
        wrong_mode = "bananas"

        self.assertRaises(InvalidArgumentError, validate_alignment_mode, wrong_mode)


class TestNetworkingValidation(TestCase):
    """Contains tests on the networking validation functions"""

    def test_validate_gcf_cutoffs(self):
        """Tests whether validate_gcf_cutoffs raises an exception when one of the
        cutoffs passed is lower than 0
        """
        cutoffs = [0.3, 0.5, -1.0]

        self.assertRaises(InvalidArgumentError, validate_gcf_cutoffs, cutoffs)


class TestDiagnosticsValidation(TestCase):
    """Contains tests on the diagnostics validation functions"""

    pass


class TestOutputValidation(TestCase):
    """Contains tests on the output validation functions"""

    def test_validate_output_dir_parent_not_exists(self):
        """Tests whether validate_output_dir raises an exception when a path is given
        to a directory of which the parent does not exist
        """
        orphan_path = Path("test/test_data/nonexistent_directory/orphan_path")

        self.assertRaises(InvalidArgumentError, validate_output_dir, orphan_path)

    def test_validate_output_dir_mkdir(self):
        """Tests whether validate_output_dir creates the output directory if necessary"""

        random_dir_str = str(randint(100000, 999999))
        random_tmp_path = Path("test/test_data/tmp") / Path(random_dir_str)

        validate_output_dir(random_tmp_path)

        self.assertTrue(random_tmp_path.exists())

        random_tmp_path.rmdir()

    def test_validate_db_return_valid(self):
        """Tests whether validate_db_path correctly returns the given db_path if it is
        not None and otherwise passes validation
        """
        # make sure tmp dir exists
        Path("test/test_data/tmp").mkdir(exist_ok=True)

        db_path = Path("test/test_data/tmp/data_sqlite.db")

        default_path = Path("test/test_data/tmp/default_db_path")

        expected_result = Path("test/test_data/tmp/data_sqlite.db")

        actual_result = validate_db_path(db_path, default_path)

        self.assertEqual(expected_result, actual_result)

    def test_validate_db_path_return_default(self):
        """Tests whether validate_db_path correctly returns the default path if both
        db_path is set to none
        """
        # make sure tmp dir exists
        Path("test/test_data/tmp").mkdir(exist_ok=True)

        db_path = None

        default_path = Path("test/test_data/tmp/default_db_path")

        expected_result = Path("test/test_data/tmp/default_db_path/data_sqlite.db")

        actual_result = validate_db_path(db_path, default_path)

        self.assertEqual(expected_result, actual_result)

    def test_validate_db_path_parent_not_exists(self):
        """Tests whether validate_db_path raises an InvalidArgumentError if a path is
        specified that points to a location without a parent directory
        """
        db_path = Path("test/test_data/nonexistent_path/data_sqlite.db")

        default_path = Path("test/test_data/tmp/default_db_path")

        self.assertRaises(InvalidArgumentError, validate_db_path, db_path, default_path)

    def test_validate_db_path_is_folder(self):
        """Tests whether validate_db_path raises an InvalidArgumentError if a path is
        specified that points to a location without a parent directory
        """
        db_path = Path("test/test_data/database")

        default_path = Path("test/test_data/tmp/default_db_path")

        self.assertRaises(InvalidArgumentError, validate_db_path, db_path, default_path)

    def test_validate_log_return_valid(self):
        """Tests whether validate_db_path correctly returns the given db_path if it is
        not None and otherwise passes validation
        """
        # make sure tmp dir exists
        Path("test/test_data/tmp").mkdir(exist_ok=True)

        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        log_path = Path("test/test_data/tmp/test.log")

        default_path = Path("test/test_data/tmp/default_log_path")

        expected_result = Path("test/test_data/tmp/test.log")

        actual_result = validate_log_path(log_path, default_path, timestamp)

        self.assertEqual(expected_result, actual_result)

    def test_validate_log_path_return_default(self):
        """Tests whether validate_db_path correctly returns the default path if both
        db_path is set to none
        """
        # make sure tmp dir exists
        Path("test/test_data/tmp").mkdir(exist_ok=True)

        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        log_path = None

        default_path = Path("test/test_data/tmp/default_log_path")

        default_log = validate_log_path(log_path, default_path, timestamp)

        # this regular expression matches the default log path + a timestamp in the
        # format: YYYY-MM-DD-HH-MM-SS + ".log"
        regexp = r"test/test_data/tmp/default_log_path/[\d-]{10}_[\d-]{8}\.log"

        self.assertRegex(str(default_log), regexp)

    def test_validate_log_path_parent_not_exists(self):
        """Tests whether validate_db_path raises an InvalidArgumentError if a path is
        specified that points to a location without a parent directory
        """

        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        log_path = Path("test/test_data/nonexistent_path/test.log")

        default_path = Path("test/test_data/tmp/default_log_path")

        self.assertRaises(
            InvalidArgumentError, validate_log_path, log_path, default_path, timestamp
        )

    def test_validate_log_path_is_folder(self):
        """Tests whether validate_db_path raises an InvalidArgumentError if a path is
        specified that points to a location without a parent directory
        """
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        db_path = Path("test/test_data/database")

        default_path = Path("test/test_data/tmp/default_log_path")

        self.assertRaises(
            InvalidArgumentError, validate_log_path, db_path, default_path, timestamp
        )
