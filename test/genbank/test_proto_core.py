"""Contains tests for the GBK class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import ProtoCore
from src.errors import InvalidGBKError
from src.data import DB


class TestProtocore(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_create_proto_core(self):
        """Tests whether a proto_core is instantiated correctly"""

        expected_number = 1

        proto_core = ProtoCore(expected_number)

        self.assertIsInstance(proto_core, ProtoCore)

    def test_parse_number(self):
        """Tests whether a protocluster/proto_core number is correctly parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")

        expected_number = 1

        feature.qualifiers["protocluster_number"] = [str(expected_number)]

        proto_core = ProtoCore.parse(feature)

        self.assertEqual(expected_number, proto_core.number)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a protocluster_number qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")

        self.assertRaises(InvalidGBKError, ProtoCore.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(InvalidGBKError, ProtoCore.parse, feature)

    def test_save(self):
        """Tests whether a GBK object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")

        expected_number = 1

        feature.qualifiers["protocluster_number"] = [str(expected_number)]

        proto_core = ProtoCore.parse(feature)

        proto_core.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_save_all(self):
        pass
