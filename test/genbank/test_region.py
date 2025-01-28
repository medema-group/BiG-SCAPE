"""Contains tests for the GBK class and functions"""

# from python
from unittest import TestCase
from pathlib import Path

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from big_scape.genbank import Region, CandidateCluster, ProtoCluster, ProtoCore, GBK
from big_scape.errors import InvalidGBKError
from big_scape.data import DB
from big_scape.enums import SOURCE_TYPE
import big_scape.enums as bs_enums


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

        region = Region.parse_as5(feature)

        self.assertEqual(expected_number, region.number)

    def test_parse_cluster_number(self):
        """Tests wether an as4 cluster number is correclty parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="cluster")
        feature.qualifiers["product"] = ["NRPS"]

        feature.qualifiers["note"] = ["Cluster number: 1"]
        expected_number = 1

        region = Region.parse_as4(feature)

        self.assertEqual(expected_number, region.number)

    def test_parse_cluster_no_number(self):
        """Tests whether parse correclty throwns an error when given a feature lacking a cluster number"""
        feature = SeqFeature(FeatureLocation(0, 100), type="cluster")
        feature.qualifiers["product"] = ["NRPS"]
        feature.qualifiers["note"] = ["Another note"]

        self.assertRaises(InvalidGBKError, Region.parse_as4, feature)

    def test_parse_region_no_number(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a region_number qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["candidate_cluster_numbers"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        self.assertRaises(InvalidGBKError, Region.parse_as5, feature)

    def test_parse_no_cand_clusters(self):
        """Tests whether parse correctly throws an error when given a feature
        lacking a candidate_cluster_numbers qualifier
        """
        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["region_number"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        gbk_file_path = Path(
            "test/test_data/metagenome_valid_gbk_input/as5_metagenome_valid...region001.gbk"
        )
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "include_categories": set(),
            "exclude_categories": set(),
        }

        parent_gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertRaises(InvalidGBKError, Region.parse_as5, feature, parent_gbk)

    def test_parse_as5_no_product(self):
        """Tests whether parse correctly throws an error when given a feature lacking a product qualifier"""

        feature = SeqFeature(FeatureLocation(0, 100), type="region")
        feature.qualifiers["region_number"] = ["1"]

        self.assertRaises(InvalidGBKError, Region.parse_as5, feature)

    def test_parse_candidate_cluster_numbers(self):
        """Tests whether a region cand_clusters is correctly parsed from a feature"""
        feature = SeqFeature(FeatureLocation(0, 100), type="region")

        expected_number = 1
        expected_cand_clusters = {1: None}

        feature.qualifiers["region_number"] = [str(expected_number)]
        feature.qualifiers["candidate_cluster_numbers"] = ["1"]
        feature.qualifiers["product"] = ["NRPS"]

        region = Region.parse_as5(feature)

        self.assertEqual(expected_cand_clusters, region.cand_clusters)

    def test_parse_as5_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a
        feature of a wrong type
        """

        feature = SeqFeature(FeatureLocation(0, 100), type="CDS")

        self.assertRaises(InvalidGBKError, Region.parse_as5, feature)

    def test_parse_as4_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a
        feature of a wrong type
        """

        feature = SeqFeature(FeatureLocation(0, 100), type="CDS")

        self.assertRaises(InvalidGBKError, Region.parse_as4, feature)

    def test_add_candidate_cluster(self):
        """Tests whether a candidate cluster is correctly added to this region"""
        region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
        region_feature.qualifiers = {
            "region_number": ["1"],
            "candidate_cluster_numbers": ["1"],
            "product": ["NRPS"],
        }

        region = Region.parse_as5(region_feature)

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

        self.assertEqual(candidate_cluster, region.cand_clusters[1])

    def test_save(self):
        """Tests whether a Region object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
        region_feature.qualifiers = {
            "region_number": ["1"],
            "candidate_cluster_numbers": ["1"],
            "product": ["NRPS"],
        }

        region = Region.parse_as5(region_feature)

        region.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_set_record_category(self):
        """Tests whether region category is set correctly"""
        core = ProtoCore(None, 1, 0, 100, False, "T1PKS", "PKS")
        pclust1 = ProtoCluster(None, 1, 0, 100, False, "T1PKS", {1: core}, "PKS")
        pclust2 = ProtoCluster(None, 2, 200, 300, False, "NRPS", {}, "NRPS")
        cand_cluster = CandidateCluster(
            None, 1, 0, 100, False, "T1PKS.NRPS", "", {1: pclust1, 2: pclust2}
        )
        # set cand_cluster category, based on protocluster categories
        cand_cluster.set_record_category(cand_cluster.get_categories())
        self.assertEqual(cand_cluster.get_category(), "NRPS.PKS")

        region = Region(None, 1, 0, 1000, False, "T1PKS.NRPS")
        region.cand_clusters = {1: cand_cluster}

        # region category is None, but can be fetched from cand_cluster
        cats = region.get_categories()
        self.assertIsNone(region.category)
        self.assertEqual(set(["PKS", "NRPS"]), cats)

        # after setting region category is equal to its cand_cluster
        region.set_record_category(cats)
        self.assertEqual(region.get_category(), "NRPS.PKS")

    def test_get_category(self):
        """Tests whether a Region category is correctly fetched"""
        region_nrps = Region(None, 1, 0, 100, None, "NRPS", "NRPS")

        self.assertEqual("NRPS", region_nrps.get_category())

        region_none = Region(None, 1, 0, 100, None, "other", None)

        self.assertEqual("Categoryless", region_none.get_category())

        region_categoryless = Region(None, 1, 0, 10, None, "", "Categoryless")

        self.assertEqual("Categoryless", region_categoryless.get_category())

    def test_as5_set_and_get_category(self):
        """Tests whether category is correctly set and read in as5 region"""
        # AS5 does not have categories
        core = ProtoCore(None, 1, 0, 100, False, "T1PKS")
        pclust1 = ProtoCluster(None, 1, 0, 100, False, "T1PKS", {1: core})
        pclust2 = ProtoCluster(None, 2, 200, 300, False, "NRPS", {})
        cand_cluster = CandidateCluster(
            None, 1, 0, 100, False, "T1PKS.NRPS", "", {1: pclust1, 2: pclust2}
        )
        region = Region(None, 1, 0, 1000, False, "T1PKS.NRPS")
        region.cand_clusters = {1: cand_cluster}

        # no categories in the region
        cats = region.get_categories()
        self.assertEqual(set([]), cats)

        # after setting, category should remain None
        region.set_record_category(cats)
        self.assertIsNone(region.category)
        self.assertEqual(region.get_category(), "Categoryless")

    def test_as4_set_and_get_category(self):
        """Tests whether category is correctly read in as4 region"""
        # AS4 does not have categories and there are no protoclusters to access
        as4_region = Region(None, 1, 0, 100, None, "other")

        self.assertEqual(as4_region.get_categories(), set([]))
        self.assertIsNone(as4_region.category)
        self.assertEqual(as4_region.get_category(), "Categoryless")
