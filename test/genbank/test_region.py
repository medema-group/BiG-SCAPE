"""Contains tests for the GBK class and functions"""
# from python
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import Region, CandidateCluster
from src.errors import InvalidGBKError
from src.data import DB


class TestRegion(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_parse_number(self):
        """Tests whether a region number is correctly parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["product"] = ["NRPS"]

        expected_number = 1

        feature.qualifiers["region_number"] = [str(expected_number)]
        feature.qualifiers["candidate_cluster_numbers"] = ["1"]

        region = Region.parse(feature)

        self.assertEqual(expected_number, region.number)

    def test_parse_cluster_number(self):
        """Tests wether an as4 cluster number is correclty parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="cluster")
        feature.qualifiers["product"] = ["NRPS"]

        feature.qualifiers["note"] = ["Cluster number: 1"]
        expected_number = 1

        region = Region.parse(feature)

        self.assertEqual(expected_number, region.number)

    def test_parse_cluster_no_number(self):
        """Tests whether parse correclty throwns an error when given a feature lacking a cluster number"""
        feature = SeqFeature(FeatureLocation(0, 100), type="cluster")
        feature.qualifiers["product"] = ["NRPS"]
        feature.qualifiers["note"] = ["Another note"]

        self.assertRaises(InvalidGBKError, Region.parse, feature)

    def test_parse_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a region_number qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["candidate_cluster_numbers"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        self.assertRaises(InvalidGBKError, Region.parse, feature)

    def test_parse_no_cand_clusters(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a candidate_cluster_numbers qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["region_number"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        self.assertRaises(InvalidGBKError, Region.parse, feature)

    def test_parse_no_product(self):
        """Tests whether parse correctly throws an error when given a feature lacking a product qualifier"""

        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["region_number"] = ["1"]

        self.assertRaises(InvalidGBKError, Region.parse, feature)

    def test_parse_candidate_cluster_numbers(self):
        """Tests whether a region cand_clusters is correctly parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="region")

        expected_number = 1
        expected_cand_clusters = {1: None}

        feature.qualifiers["region_number"] = [str(expected_number)]
        feature.qualifiers["candidate_cluster_numbers"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        region = Region.parse(feature)

        self.assertEqual(expected_cand_clusters, region.cand_clusters)

    def test_parse_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a
        feature of a wrong type
        """

        feature = SeqFeature(FeatureLocation(0, 100), type="CDS")

        self.assertRaises(InvalidGBKError, Region.parse, feature)

    def test_add_candidate_cluster(self):
        """Tests whether a candidate cluster is correctly added to this region"""
        region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
        region_feature.qualifiers = {
            "region_number": ["1"],
            "candidate_cluster_numbers": ["1"],
            "product": ["NRPS"],
        }

        region = Region.parse(region_feature)

        candidate_cluster_feature = SeqFeature(
            FeatureLocation(0, 100), type="cand_cluster"
        )
        candidate_cluster_feature.qualifiers = {
            "candidate_cluster_number": ["1"],
            "kind": ["neighbouring"],
            "protoclusters": ["1"],
            "product": ["NRPS"],
        }

        candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

        region.add_cand_cluster(candidate_cluster)

    def test_save(self):
        """Tests whether a Region object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
        region_feature.qualifiers = {
            "region_number": ["1"],
            "candidate_cluster_numbers": ["1"],
            "product": ["NRPS"],
        }

        region = Region.parse(region_feature)

        region.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
