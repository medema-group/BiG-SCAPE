"""Contains tests for the Protocluster class and functions"""

from unittest import TestCase

from Bio.SeqFeature import SeqFeature

from src.genbank.protocluster import Protocluster
from src.genbank.proto_core import Protocore
from src.errors.genbank import InvalidGBKError


class TestProtocluster(TestCase):
    """Test class for GBK Protocluster parsing tests"""

    def test_create_protocluster(self):
        """Tests whether a Protocluster is instatiated correclty"""

        expected_number = 1

        protocluster = Protocluster(expected_number)

        self.assertIsInstance(protocluster, Protocluster)

    def test_parse_number(self):
        """Tests whether a Protocluster number is correctly parsed from a feature"""

        expected_number = 1

        protocluster_feature = SeqFeature(type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
        }

        protocluster = Protocluster.parse(protocluster_feature)

        self.assertEqual(expected_number, protocluster.number)

    def test_parse_category(self):
        """Tests whether a Protocluster category is correctly parsed from a feature"""

        expected_category = "NRPS"

        protocluster_feature = SeqFeature(type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
        }

        protocluster = Protocluster.parse(protocluster_feature)

        self.assertEqual(expected_category, protocluster.category)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a protocluster_number qualifier
        """
        feature = SeqFeature(type="protocluster")
        expected_category = "NRPS"
        feature.qualifiers["category"] = [expected_category]

        self.assertRaises(InvalidGBKError, Protocluster.parse, feature)

    def test_parse_no_category(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a category qualifier
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

    def test_add_protocore(self):
        """Tests whether a protocore is correctly added to this protocluster"""

        protocluster_feature = SeqFeature(type="protocluster")
        protocluster_feature.qualifiers = {
            "protocluster_number": ["1"],
            "category": ["NRPS"],
        }

        protocluster = Protocluster.parse(protocluster_feature)

        protocore_feature = SeqFeature(type="proto_core")
        protocore_feature.qualifiers = {
            "protocluster_number": ["1"],
        }

        protocore = Protocore.parse(protocore_feature)

        protocluster.add_proto_core(protocore)
