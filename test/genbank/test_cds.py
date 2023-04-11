"""Contains tests for the CDS class and functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import GBK, CDS
from src.errors import InvalidGBKError
from src.data import DB


class TestCDS(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_create_cds(self):
        """Tests whether a cds is instatiated correclty"""

        expected_start = 1
        expected_end = 5

        cds = CDS(expected_start, expected_end)

        self.assertIsInstance(cds, CDS)

    def test_parse_aa_seq(self):
        """Tests whether an aa sequence is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["gene_kind"] = "biosynthetic"

        expected_aa_seq = "MDRAAGMSVRLK"

        feature.qualifiers["translation"] = [expected_aa_seq]

        cds = CDS.parse(feature)

        self.assertEqual(expected_aa_seq, cds.aa_seq)

    def test_parse_no_aa_seq(self):
        """Tests whehter parse correctly throws an error when given a feature
        lacking a translation qualifier
        """

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["gene_kind"] = ["biosynthetic"]

        self.assertRaises(InvalidGBKError, CDS.parse, feature)

    def test_parse_gene_kind(self):
        """Tests whether a gene_kind is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["translation"] = ["MDRAAGMSVRLK"]

        expected_gene_kind = "biosynthetic"

        feature.qualifiers["gene_kind"] = [expected_gene_kind]

        cds = CDS.parse(feature)

        self.assertEqual(expected_gene_kind, cds.gene_kind)

    def test_parse_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="region")

        self.assertRaises(InvalidGBKError, CDS.parse, feature)

    def test_save_cds(self):
        """Tests whether a CDS can be correctly saved to the database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, "query")

        gbk.save_all()

        DB.commit()

        DB.save_to_disk(Path("tmp/db.db"))

        # 1 gbk, 11 bgc records
        expected_row_count = 37

        actual_row_count = 0

        cursor_result = DB.execute_raw_query("SELECT * FROM cds;")
        actual_row_count += len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
