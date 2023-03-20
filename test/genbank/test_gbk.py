"""Contains tests for the GBK class and functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.Seq import Seq

# from other modules
from src.genbank import GBK, Region, ProtoCore, CDS
from src.errors import InvalidGBKError
from src.data import DB


class TestGBK(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

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

    def test_populate_cds(self):
        """Tests whether parsing a GBK correctly populates the underlying CDSs"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        self.assertIsInstance(gbk.genes[0], CDS)

    def test_populate_hierarchical_objects(self):
        """Tests whether parsing a GBK correclty generates parent-child feature relations
        via checking for presence of the lowest level child - proto_core"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        proto_core = gbk.region.cand_clusters[1].proto_clusters[1].proto_core[1]

        self.assertIsInstance(proto_core, ProtoCore)

    def test_parse_gbk_has_dna_seq(self):
        """Tests whether parsing a GBK correclty has DNA sequence"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        dna_sequence = gbk.nt_seq

        self.assertIsInstance(dna_sequence, Seq)

    def test_save(self):
        """Tests whether a GBK object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk = GBK.parse(gbk_file_path)

        gbk.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM gbk;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_save_all(self):
        pass
