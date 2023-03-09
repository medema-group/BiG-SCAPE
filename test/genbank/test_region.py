"""Contains tests for the GBK class and functions"""
from unittest import TestCase

from Bio.SeqFeature import SeqFeature

from src.genbank.region import Region


class TestRegion(TestCase):
    """Test class for base GBK parsing tests"""

    def test_create_region(self):
        """Tests whether a region is instantiated correctly"""
        feature = SeqFeature(type="region")
        feature.qualifiers["region_number"] = ["1"]

        region = Region.create_region(feature)

        self.assertIsInstance(region, Region)

    def test_create_region_number(self):
        """Tests whether a region is instantiated correctly"""
        feature = SeqFeature(type="region")

        expected_number = 1

        feature.qualifiers["region_number"] = [str(expected_number)]

        region = Region.create_region(feature)

        self.assertEqual(expected_number, region.number)

    def test_create_region_no_number(self):
        """Tests whether create_region correctly throws an error when given a feature
        lacking a region_number qualifier
        """
        feature = SeqFeature(type="region")

        self.assertRaises(ValueError, Region.create_region, feature)

    def test_create_region_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="CDS")

        self.assertRaises(ValueError, Region.create_region, feature)
