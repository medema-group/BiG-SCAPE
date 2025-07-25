"""module containing unit tests for loading input files"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
import big_scape.enums as bs_enums
from big_scape.dereplicating.gbk_components.gbk import GBK
from big_scape.dereplicating.gbk_components.cds import CDS
from big_scape.errors import InvalidGBKError

# from this module
from big_scape.dereplicating.input_data_loading import (
    load_input_folder,
    parse_gbk_files,
    gbk_factory,
    get_parser_functions,
    parse_seqIO,
    load_input_data,
    make_gbk_name
)


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

        gbk_data = parse_gbk_files(input_gbk_paths, bs_enums.SOURCE_TYPE.QUERY, run)

        gbk_data_list = list(gbk_data)

        self.assertIsInstance(gbk_data_list[0], tuple)
        self.assertIsInstance(gbk_data_list[0][0], str)
        self.assertIsInstance(gbk_data_list[0][1], Path)
        self.assertIsInstance(gbk_data_list[0][2], str)
        self.assertIsInstance(gbk_data_list[0][3], SeqRecord)

    def test_make_gbk_name(self):
        """tests whether the make_gbk_name function correctly creates a GBK name"""

        run = {
            "input_dir": Path("test/test_data/"),
        }

        gbk_path_abs = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk_hash = "hash"

        gbk_name = make_gbk_name(run, gbk_path_abs, gbk_hash)

        expected_name = "valid_gbk_folder.valid_input_region.gbk.hash"

        self.assertIsInstance(gbk_name, str)
        self.assertEqual(gbk_name, expected_name)

    def test_get_parser_functions(self):
        """Tests whether the get_parser_functions function correctly returns the correct parser functions"""

        run = {
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
        }

        parser_functions = get_parser_functions(run)

        self.assertIsInstance(parser_functions, dict)
        self.assertIsInstance(
            parser_functions[bs_enums.FEATURE_TYPE.CDS], type(CDS.parse)
        )

    def test_parse_seqIO_record(self):
        """Tests whether the gbk_factory function correctly creates a GBK object"""

        gbk = GBK("name", Path("test_path"), "hash", 10, "1", bs_enums.SOURCE_TYPE.QUERY)

        nt_seq = Seq("ATGCAGCAGGACGGCACACAGCAGGACCGGATCAAGCAGAGTCCCGCCCCTCTCTGA")
        seqIO_record = SeqRecord(id="test", seq=nt_seq)

        expected_transl_nt_seq = Seq("MQQDGTQQDRIKQSPAPL")

        feature = SeqFeature(
            FeatureLocation(5, 10, strand=1),
            type="CDS",
            qualifiers={"translation": expected_transl_nt_seq},
        )

        seqIO_record.features = [feature]

        gbk = parse_seqIO(gbk, seqIO_record, bs_enums.RUN_MODE.DEREPLICATE)

        self.assertIsInstance(gbk, GBK)
        self.assertEqual(len(gbk.components), 1)

        components_types = list(gbk.components.keys())
        self.assertEqual(len(set(components_types)), 1)
        self.assertEqual(components_types[0], bs_enums.COMPONENTS.CDS)

        components_values = list(gbk.components.values())
        self.assertEqual(len(components_values), 1)
        self.assertIsInstance(components_values[0][0], CDS)

    def test_gbk_factory_dereplicate_mode(self):
        """Tests whether the gbk_factory function correctly creates a GBK object in dereplicate mode"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
        }

        input_gbk_paths = load_input_folder(run)

        gbk_data = parse_gbk_files(input_gbk_paths, bs_enums.SOURCE_TYPE.QUERY, run)

        gbk_1 = next(gbk_data)

        gbk = gbk_factory(gbk_1, run)

        # test whether only CDS components are present
        components_types = list(gbk.components.keys())
        self.assertEqual(len(set(type(x) for x in components_types)), 1)
        self.assertEqual(components_types[0], bs_enums.COMPONENTS.CDS)

    def test_gbk_factory_no_cds_error(self):
        """Tests whether the gbk_factory function correctly raises an error and returns None if no CDS present"""

        run = {
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
        }

        seqIO_record = SeqRecord(id="test", seq=Seq("ATG"))

        gbk_data = ("name", Path("path"), "hash", seqIO_record, bs_enums.SOURCE_TYPE.QUERY)

        gbk_list = gbk_factory(gbk_data, run)

        self.assertEqual(gbk_list, None)

    def test_load_input_data(self):
        """Tests whether the load_input_data function correctly loads and filters files"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
        }

        gbk_list = load_input_folder(run)

        self.assertEqual(len(gbk_list), 2)

    def test_load_input_data_empty_gbk_folder(self):
        """Tests whether the load_input_data function correctly raises an error if no GBK files are found"""

        run = {
            "input_dir": Path("test/test_data/empty_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
        }

        with self.assertRaises(FileNotFoundError):
            load_input_data(run)

    def test_load_input_data_invalid_gbk_folder(self):
        """Tests whether the load_input_data function correctly raises an error if the input path is not a directory"""

        run = {
            "input_dir": Path("test/test_data/invalid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
            "cores": 1
        }

        with self.assertRaises(InvalidGBKError):
            load_input_data(run)
