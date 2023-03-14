"""Contains tests for the Protocluster class and functions"""

from unittest import TestCase

from Bio.SeqFeature import SeqFeature

from src.genbank.protocluster import Protocluster
from src.errors.genbank import InvalidGBKError


class TestProtocluster(TestCase):
    """Test class for GBK Protocluster parsing tests"""

    def test_create_protocluster(self):
        """Tests whether a Protocluster is instatiated correclty"""

        expected_number = 1
        expected_category = "any category"
        protocluster = Protocluster(expected_number, expected_category)

        self.assertIsInstance(protocluster, Protocluster)

    def test_parse_number(self):
        """Tests whether a Protocluster number is correctly parsed from a feature"""
        feature = SeqFeature(type="protocluster")

        expected_number = 1
        expected_category = "any category"
        feature.qualifiers["protocluster_number"] = [str(expected_number)]
        feature.qualifiers["category"] = [expected_category]

        protocluster = Protocluster.parse(feature)

        self.assertEqual(expected_number, protocluster.number)

    def test_parse_category(self):
        """Tests whether a Protocluster category is correctly parsed from a feature"""
        feature = SeqFeature(type="protocluster")

        expected_number = 1
        expected_category = "any category"
        feature.qualifiers["protocluster_number"] = [str(expected_number)]
        feature.qualifiers["category"] = [expected_category]

        protocluster = Protocluster.parse(feature)

        self.assertEqual(expected_category, protocluster.category)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a protocluster_number qualifier
        """
        feature = SeqFeature(type="protocluster")
        expected_category = "any category"
        feature.qualifiers["category"] = [expected_category]

        self.assertRaises(InvalidGBKError, Protocluster.parse, feature)

    def test_parse_no_category(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a kind qualifier
        """
        feature = SeqFeature(type="protocluster")
        expected_number = 1
        feature.qualifiers["protocluster_number"] = [str(expected_number)]

        self.assertRaises(InvalidGBKError, Protocluster.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(InvalidGBKError, Protocluster.parse, feature)
