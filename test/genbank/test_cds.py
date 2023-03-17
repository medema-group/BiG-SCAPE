"""Contains tests for the CDS class and functions"""

# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import CDS
from src.errors import InvalidGBKError


class TestCDS(TestCase):
    """Test class for base GBK parsing tests"""

    def test_create_cds(self):
        """Tests whether a cds is instatiated correclty"""

        expected_start = 1
        expected_end = 5

        cds = CDS(expected_start, expected_end)

        self.assertIsInstance(cds, CDS)

    def test_parse_aa_seq(self):
        """Tests whether an aa sequence is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10), strand=1, type="CDS")
        feature.qualifiers["gene_kind"] = "biosynthetic"

        expected_aa_seq = "MDRAAGMSVRLK"

        feature.qualifiers["translation"] = [expected_aa_seq]

        cds = CDS.parse(feature)

        self.assertEqual(expected_aa_seq, cds.aa_seq)

    def test_parse_no_aa_seq(self):
        """Tests whehter parse correctly throws an error when given a feature
        lacking a translation qualifier
        """

        feature = SeqFeature(FeatureLocation(5, 10), strand=1, type="CDS")
        feature.qualifiers["gene_kind"] = ["biosynthetic"]

        self.assertRaises(InvalidGBKError, CDS.parse, feature)

    def test_parse_gene_kind(self):
        """Tests whether a gene_kind is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10), strand=1, type="CDS")
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
