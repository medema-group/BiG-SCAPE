"""Contains tests for the Protocluster class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from big_scape.genbank import (
    ProtoCluster,
    ProtoCore,
    MergedProtoCluster,
    CandidateCluster,
)
from big_scape.errors import InvalidGBKError
from big_scape.data import DB


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
        """Tests whether a protocore is correctly added to this protocluster and its category updated"""

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

        self.assertEqual(protocore.category, protocluster.category)

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

        proto_cluster.save(0)

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_merge(self):
        """tests whether two proto_clusters are correclty merged"""

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        protocluster_feature_a = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_a.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster_a = ProtoCluster.parse(protocluster_feature_a)
        protocluster_a.add_proto_core(proto_core_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocluster_feature_b = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_b.qualifiers = {
            "protocluster_number": ["2"],
            "category": ["PKS"],
            "product": ["PKS"],
        }

        protocluster_b = ProtoCluster.parse(protocluster_feature_b)
        protocluster_b.add_proto_core(proto_core_b)

        merged_protoclusters = MergedProtoCluster.merge(
            [protocluster_a, protocluster_b]
        )

        expected_number = "1_2"

        self.assertEqual(expected_number, merged_protoclusters.merged_number)

    def test_save_merged(self):
        """Tests whether a merged protocluster is saved correclty to the database"""

        DB.create_in_mem()

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        protocluster_feature_a = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_a.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster_a = ProtoCluster.parse(protocluster_feature_a)
        protocluster_a.add_proto_core(proto_core_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocluster_feature_b = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_b.qualifiers = {
            "protocluster_number": ["2"],
            "category": ["PKS"],
            "product": ["PKS"],
        }

        protocluster_b = ProtoCluster.parse(protocluster_feature_b)
        protocluster_b.add_proto_core(proto_core_b)

        merged_protoclusters = MergedProtoCluster.merge(
            [protocluster_a, protocluster_b]
        )
        merged_protoclusters.save(0)

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_load_merged(self):
        """tests whether a merged protocluster is loaded correctly from the database"""

        DB.create_in_mem()

        feature_a = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_a.qualifiers["protocluster_number"] = [str(1)]
        feature_a.qualifiers["product"] = ["NRPS"]
        proto_core_a = ProtoCore.parse(feature_a)

        protocluster_feature_a = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_a.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        protocluster_a = ProtoCluster.parse(protocluster_feature_a)
        protocluster_a.add_proto_core(proto_core_a)

        feature_b = SeqFeature(FeatureLocation(0, 100), type="proto_core")
        feature_b.qualifiers["protocluster_number"] = [str(2)]
        feature_b.qualifiers["product"] = ["PKS"]
        proto_core_b = ProtoCore.parse(feature_b)

        protocluster_feature_b = SeqFeature(
            FeatureLocation(0, 100), type="protocluster"
        )
        protocluster_feature_b.qualifiers = {
            "protocluster_number": ["2"],
            "category": ["PKS"],
            "product": ["PKS"],
        }

        protocluster_b = ProtoCluster.parse(protocluster_feature_b)
        protocluster_b.add_proto_core(proto_core_b)

        merged_protoclusters = MergedProtoCluster.merge(
            [protocluster_a, protocluster_b]
        )
        merged_protoclusters.save(0)

        candidate_cluster_feature = SeqFeature(
            FeatureLocation(0, 100), type="cand_cluster"
        )
        candidate_cluster_feature.qualifiers = {
            "candidate_cluster_number": ["1"],
            "kind": ["single"],
            "protoclusters": ["1"],
            "product": ["NRPS"],
        }

        cand_cluster = CandidateCluster.parse(candidate_cluster_feature)
        cand_cluster_dict = {0: cand_cluster}

        ProtoCluster.load_all(cand_cluster_dict)

        loaded_proto_cluster = cand_cluster.proto_clusters["1_2"]

        self.assertIsInstance(loaded_proto_cluster, MergedProtoCluster)
