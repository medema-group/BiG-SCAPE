"""Contains tests for the Protocluster class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import ProtoCluster, ProtoCore
from src.errors import InvalidGBKError
from src.data import DB


class TestProtocluster(TestCase):
    """Test class for GBK Protocluster parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_parse_number(self):
        """Tests whether a Protocluster number is correctly parsed from a feature"""

        expected_number = 1

        protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster = ProtoCluster.parse(protocluster_feature)

        self.assertEqual(expected_number, protocluster.number)

    def test_parse_category(self):
        """Tests whether a Protocluster category is correctly parsed from a feature"""

        expected_category = "NRPS"

        protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster = ProtoCluster.parse(protocluster_feature)

        self.assertEqual(expected_category, protocluster.category)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a protocluster_number qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        expected_category = "NRPS"
        feature.qualifiers["category"] = [expected_category]

        self.assertRaises(InvalidGBKError, ProtoCluster.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(FeatureLocation(0, 100), type="CDS")

        self.assertRaises(InvalidGBKError, ProtoCluster.parse, feature)

    def test_add_protocore(self):
        """Tests whether a protocore is correctly added to this protocluster"""

        protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster = ProtoCluster.parse(protocluster_feature)

        protocore_feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        protocore_feature.qualifiers = {
            "protocluster_number": ["1"],
            "product": ["NRPS"],
        }

        protocore = ProtoCore.parse(protocore_feature)

        protocluster.add_proto_core(protocore)

    def test_save(self):
        """Tests whether a ProtoCluster object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        proto_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        proto_cluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        proto_cluster = ProtoCluster.parse(proto_cluster_feature)

        proto_cluster.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
