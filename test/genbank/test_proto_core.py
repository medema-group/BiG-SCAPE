"""Contains tests for the GBK class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from big_scape.genbank import ProtoCore, MergedProtoCore, ProtoCluster
from big_scape.errors import InvalidGBKError
from big_scape.data import DB


class TestProtocore(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_parse_number(self):
        """Tests whether a protocluster/proto_core number is correctly parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")

        expected_number = 1

        feature.qualifiers["protocluster_number"] = [str(expected_number)]
        feature.qualifiers["product"] = ["NRPS"]

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
        feature.qualifiers["product"] = ["NRPS"]

        proto_core = ProtoCore.parse(feature)

        proto_core.save(0)

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_merge(self):
        """Tests whether two protocores are correclty merged"""

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocores = [proto_core_a, proto_core_b]
        merged = MergedProtoCore.merge(protocores)

        expected_number = "1_2"

        self.assertEqual(expected_number, merged.merged_number)

    def test_save_merged(self):
        """Tests whether a merged protocore is saved correclty to the database"""

        DB.create_in_mem()

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocores = [proto_core_a, proto_core_b]
        merged = MergedProtoCore.merge(protocores)

        merged.save(0)

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_load_merged(self):
        """Tests whether a merged protocore is loaded correctly from the database"""

        DB.create_in_mem()

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocores = [proto_core_a, proto_core_b]
        merged = MergedProtoCore.merge(protocores)

        merged.save(0)

        protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["0"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster = ProtoCluster.parse(protocluster_feature)
        protocluster_dict = {0: protocluster}

        ProtoCore.load_all(protocluster_dict)

        loaded_protocore = protocluster.proto_core["1_2"]

        self.assertIsInstance(loaded_protocore, MergedProtoCore)
