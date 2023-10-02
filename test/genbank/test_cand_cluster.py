"""Contains tests for the Cand_cluster class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation


# from other modules
from big_scape.genbank import CandidateCluster, ProtoCluster
from big_scape.errors import InvalidGBKError
from big_scape.data import DB


class TestCandidateCluster(TestCase):
    """Test class for GBK Cand_cluster parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_parse_number(self):
        """Tests whether a CandidateCluster number is correctly parsed from a feature"""

        expected_number = 1
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

        self.assertEqual(expected_number, cand_cluster.number)

    def test_parse_kind(self):
        """Tests whether a CandidateCluster kind is correctly parsed from a feature"""

        expected_kind = "single"
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

        self.assertEqual(expected_kind, cand_cluster.kind)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a cand_cluster_number qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)

    def test_parse_no_kind(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a kind qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
        expected_number = 1
        feature.qualifiers["candidate_cluster_number"] = [str(expected_number)]

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)

    def test_add_protocluster(self):
        """Tests whether a protocluster is correctly added to this candidate cluster"""

        candidate_cluster_feature = SeqFeature(
            FeatureLocation(0, 100), type="cand_cluster"
        )
        candidate_cluster_feature.qualifiers = {
            "candidate_cluster_number": ["1"],
            "kind": ["single"],
            "protoclusters": ["1"],
            "product": ["NRPS"],
        }

        candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

        proto_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
        proto_cluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
            "product": ["NRPS"],
        }

        proto_cluster = ProtoCluster.parse(proto_cluster_feature)

        candidate_cluster.add_proto_cluster(proto_cluster)

    def test_save(self):
        """Tests whether a CandidateCluster object is correctly stored in the SQLite
        database
        """

        DB.create_in_mem()

        candidate_cluster_feature = SeqFeature(
            FeatureLocation(0, 100), type="cand_cluster"
        )
        candidate_cluster_feature.qualifiers = {
            "candidate_cluster_number": ["1"],
            "kind": ["single"],
            "protoclusters": ["1"],
            "product": ["NRPS"],
        }

        candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

        candidate_cluster.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
