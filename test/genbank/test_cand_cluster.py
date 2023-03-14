"""Contains tests for the Cand_cluster class and functions"""

from unittest import TestCase

from Bio.SeqFeature import SeqFeature

from src.genbank.cand_cluster import CandidateCluster
from src.errors.genbank import InvalidGBKError


class TestCandidateCluster(TestCase):
    """Test class for GBK Cand_cluster parsing tests"""

    def test_create_cand_cluster(self):
        """Tests whether a CandidateCluster is instatiated correclty"""

        expected_number = 1
        expected_kind = "any kind"
        cand_cluster = CandidateCluster(expected_number, expected_kind)

        self.assertIsInstance(cand_cluster, CandidateCluster)

    def test_parse_number(self):
        """Tests whether a CandidateCluster number is correctly parsed from a feature"""
        feature = SeqFeature(type="cand_cluster")

        expected_number = 1
        expected_kind = "any kind"
        feature.qualifiers["candidate_cluster_number"] = [str(expected_number)]
        feature.qualifiers["kind"] = [str(expected_kind)]

        cand_cluster = CandidateCluster.parse(feature)

        self.assertEqual(expected_number, cand_cluster.number)

    def test_parse_kind(self):
        """Tests whether a CandidateCluster kind is correctly parsed from a feature"""
        feature = SeqFeature(type="cand_cluster")

        expected_number = 1
        expected_kind = "any kind"
        feature.qualifiers["candidate_cluster_number"] = [str(expected_number)]
        feature.qualifiers["kind"] = [expected_kind]

        cand_cluster = CandidateCluster.parse(feature)

        self.assertEqual(expected_kind, cand_cluster.kind)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a cand_cluster_number qualifier
        """
        feature = SeqFeature(type="cand_cluster")
        expected_kind = "any kind"
        feature.qualifiers["kind"] = [expected_kind]

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)

    def test_parse_no_kind(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a kind qualifier
        """
        feature = SeqFeature(type="cand_cluster")
        expected_number = 1
        feature.qualifiers["candidate_cluster_number"] = [str(expected_number)]

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(InvalidGBKError, CandidateCluster.parse, feature)
