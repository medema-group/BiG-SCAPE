"""Contains tests for the GBK class and functions"""
from unittest import TestCase

from Bio.SeqFeature import SeqFeature

from src.genbank.proto_core import ProtoCore
from src.errors.genbank import InvalidGBKError


class TestProtocore(TestCase):
    """Test class for base GBK parsing tests"""

    def test_create_proto_core(self):
        """Tests whether a proto_core is instantiated correctly"""

        expected_number = 1

        proto_core = ProtoCore(expected_number)

        self.assertIsInstance(proto_core, ProtoCore)

    def test_parse_number(self):
        """Tests whether a protocluster/proto_core number is correctly parsed from a feature"""
        feature = SeqFeature(type="proto_core")

        expected_number = 1

        feature.qualifiers["protocluster_number"] = [str(expected_number)]

        proto_core = ProtoCore.parse(feature)

        self.assertEqual(expected_number, proto_core.number)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a protocluster_number qualifier
        """
        feature = SeqFeature(type="proto_core")

        self.assertRaises(InvalidGBKError, ProtoCore.parse, feature)

    def test_parse_wrong_type(self):
        """Tests whether parse correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(InvalidGBKError, ProtoCore.parse, feature)
