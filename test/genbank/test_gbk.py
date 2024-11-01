"""Contains tests for the GBK class and functions"""

# from python
import logging
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.Seq import Seq

# from other modules
from big_scape.genbank import GBK, Region, ProtoCore, CDS
from big_scape.errors import InvalidGBKError
from big_scape.data import DB
from big_scape.enums import SOURCE_TYPE
import big_scape.enums as bs_enums

import big_scape.genbank as bs_gbk


class TestGBK(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_parse_gbk(self):
        """Tests whether a GBK is instantiated correctly"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk, GBK)

    def test_parse_as4gbk(self):
        """Tests whether an as4 GBK is instantiated correctly"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/CM001015.1.cluster001.gbk"
        )

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "force_gbk": False,
            "include_categories": set(),
            "exclude_categories": set(),
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk, GBK)

    def test_parse_as4gbk_product(self):
        """Tests whether an as4 GBK's product is parsed correclty"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/CM000578.1.cluster042.gbk"
        )

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "force_gbk": False,
            "include_categories": set(),
            "exclude_categories": set(),
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)
        expected_region_product = "nrps.t1pks"
        self.assertEqual(gbk.region.product, expected_region_product)

    def test_error_parse_as4_gbk_classify_legacy_weights(self):
        """Tests whether an error is raise when trying to parse a as4 GBK with legacy weights"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/CM001015.1.cluster001.gbk"
        )
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": True,
            "legacy_weights": True,
        }
        with self.assertRaises(InvalidGBKError):
            GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

    def test_parse_as6_metagenome_gbk_classify_legacy_weights(self):
        """Tests whether an as6 metagenome GBK is parsed with no errors
        even when legacy weights are used
        """

        gbk_file_path = Path(
            "test/test_data/metagenome_valid_gbk_input/as6_metagenome_valid...region001.gbk"
        )
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.FLAT,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": True,
            "legacy_weights": True,
        }

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk, GBK)

    def test_get_sub_records(self):
        """Tests whether all records are correctly retrieved from a gbk, internally also testing that collapse
        protoclusters works correclty"""

        # a neighbouring cand_cluster, one interleaved and one hybrid cand_clusters
        gbk_file_path_1 = Path(
            "test/test_data/metagenome_valid_gbk_input/JCM_4504.region28.gbk"
        )
        # a neighbouring cand_cluster, 3 single cand_clusters
        gbk_file_path_2 = Path(
            "test/test_data/metagenome_valid_gbk_input/NS1.region05.gbk"
        )

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.FLAT,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": True,
            "legacy_weights": True,
        }

        gbk_1 = GBK.parse(gbk_file_path_1, SOURCE_TYPE.QUERY, run)
        gbk_2 = GBK.parse(gbk_file_path_2, SOURCE_TYPE.QUERY, run)

        all_protoclusters_1 = bs_gbk.bgc_record.get_sub_records(
            gbk_1.region, bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER
        )
        protocluster_numbers_1 = [pc.number for pc in all_protoclusters_1]
        protocluster_numbers_1.sort()

        all_protoclusters_2 = bs_gbk.bgc_record.get_sub_records(
            gbk_2.region, bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER
        )
        protocluster_numbers_2 = [pc.number for pc in all_protoclusters_2]
        protocluster_numbers_2.sort()

        all_proto_cores_1 = bs_gbk.bgc_record.get_sub_records(
            gbk_1.region, bs_enums.genbank.RECORD_TYPE.PROTO_CORE
        )
        proto_core_numbers_1 = [pc.number for pc in all_proto_cores_1]
        proto_core_numbers_1.sort()

        all_proto_cores_2 = bs_gbk.bgc_record.get_sub_records(
            gbk_2.region, bs_enums.genbank.RECORD_TYPE.PROTO_CORE
        )
        proto_core_numbers_2 = [pc.number for pc in all_proto_cores_2]
        proto_core_numbers_2.sort()

        all_cand_clusters_1 = bs_gbk.bgc_record.get_sub_records(
            gbk_1.region, bs_enums.genbank.RECORD_TYPE.CAND_CLUSTER
        )
        cand_cluster_numbers_1 = [cc.number for cc in all_cand_clusters_1]
        cand_cluster_numbers_1.sort()
        cc_1_1_pcs = list(gbk_1.region.cand_clusters[1].proto_clusters.keys())
        cc_1_1_pcs.sort()
        cc_1_2_pcs = list(gbk_1.region.cand_clusters[2].proto_clusters.keys())
        cc_1_2_pcs.sort()
        cc_1_3_pcs = list(gbk_1.region.cand_clusters[3].proto_clusters.keys())
        cc_1_3_pcs.sort()

        all_cand_clusters_2 = bs_gbk.bgc_record.get_sub_records(
            gbk_2.region, bs_enums.genbank.RECORD_TYPE.CAND_CLUSTER
        )
        cand_cluster_numbers_2 = [cc.number for cc in all_cand_clusters_2]
        cand_cluster_numbers_2.sort()
        cc_2_1_pcs = list(gbk_2.region.cand_clusters[1].proto_clusters.keys())
        cc_2_1_pcs.sort()
        cc_2_2_pcs = list(gbk_2.region.cand_clusters[2].proto_clusters.keys())
        cc_2_2_pcs.sort()
        cc_2_3_pcs = list(gbk_2.region.cand_clusters[3].proto_clusters.keys())
        cc_2_3_pcs.sort()
        cc_2_4_pcs = list(gbk_2.region.cand_clusters[4].proto_clusters.keys())
        cc_2_4_pcs.sort()

        seen_dict = {
            "pcl_1": protocluster_numbers_1,
            "pc_1": proto_core_numbers_1,
            "cc_1": cand_cluster_numbers_1,
            "cc_1_1": cc_1_1_pcs,
            "cc_1_2": cc_1_2_pcs,
            "cc_1_3": cc_1_3_pcs,
            "pcl_2": protocluster_numbers_2,
            "pc_2": proto_core_numbers_2,
            "cc_2": cand_cluster_numbers_2,
            "cc_2_1": cc_2_1_pcs,
            "cc_2_2": cc_2_2_pcs,
            "cc_2_3": cc_2_3_pcs,
            "cc_2_4": cc_2_4_pcs,
        }

        expected_dict = {
            "pcl_1": [1, 3],
            "pc_1": [1, 3],
            "cc_1": [1, 2, 3],
            "cc_1_1": [1, 3],
            "cc_1_2": [1],
            "cc_1_3": [3],
            "pcl_2": [1, 2, 3],
            "pc_2": [1, 2, 3],
            "cc_2": [1, 2, 3, 4],
            "cc_2_1": [1, 2, 3],
            "cc_2_2": [1],
            "cc_2_3": [2],
            "cc_2_4": [3],
        }

        self.assertEqual(seen_dict, expected_dict)

    def test_parse_as4_no_cluster_feature(self):
        """Tests whether an as4 gbk has no cluster feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as4_no_cluster_feature.gbk"
        )
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "force_gbk": False,
            "include_categories": set(),
            "exclude_categories": set(),
        }
        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY, run
        )

    def test_parse_metagenome_gbk(self):
        """Tests whether a metagenome GBK is instantiated correclty"""

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
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk, GBK)

    def test_parse_mibig_bac_region_only_gbk(self):
        """Tests whether a MIBiG bacterial gbk file, with only regions and CDSs, is instantiated correclty"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002476.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG, run)

        self.assertIsInstance(gbk, GBK)

    def test_parse_mibig_fun_region_only_gbk(self):
        """Tests whether a MIBiG fungal gbk file, with only regions and CDSs, is instantiated correclty"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002609.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG, run)

        self.assertIsInstance(gbk, GBK)

    def test_parse_gbk_multiple_regions(self):
        """Tests whether a GBK file has more than one region"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_multiple_regions_folder/valid_input_multiple_regions.gbk"
        )
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY, run
        )

    def test_parse_as5_no_region(self):
        """Tests whether a GBK file has no region feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_no_region_feature.gbk"
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
        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY, run
        )

    def test_parse_as5_no_cand_cluster(self):
        """Tests whether a GBK file with no cand_cluster features gives warning"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_no_cand_cluster_feature.gbk"
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

        with self.assertLogs(level=logging.INFO) as cm:
            logging.info("nonsense")
            GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        # cm.output a list of strings of all the logs
        str = "does not contain any cand_cluster features"
        warning = any(str in log for log in cm.output)

        self.assertEqual(warning, True)

    def test_parse_as5_missing_protocluster(self):
        """Tests whether a GBK file has missing protocluster feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_missing_protocluster_feature.gbk"
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

        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY, run
        )

    def test_parse_as5_missing_proto_core(self):
        """Tests whether a GBK file has missing proto_core feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_missing_proto_core_feature.gbk"
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

        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY, run
        )

    def test_as_invalid_version_number(self):
        """Tests whether found as version number is valid (is a number)"""
        gbk_file_path = Path("test/test_data/invalid_gbk_folder/as_invalid_version.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.MIBIG, run
        )

    def test_as_below_6_legacy_weights(self):
        """Tests valid combinations of classify and legacy_weights"""
        gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/CM001015.1.cluster001.gbk"
        )
        run = {
            "force_gbk": False,
            "classify": bs_enums.CLASSIFY_MODE.LEGACY,
            "legacy_weights": True,
            "include_categories": set(),
            "exclude_categories": set(),
        }

        # no errors for legacy + legacy_weights
        GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG, run)

        run["classify"] = bs_enums.CLASSIFY_MODE.CLASS

        # error for not legacy + legacy_weights (on AS4 gbk)
        self.assertRaises(
            InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.MIBIG, run
        )

    def test_populate_regions(self):
        """Tests whether parsing a GBK correctly populates the underlying region"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk.region, Region)

    def test_populate_cds(self):
        """Tests whether parsing a GBK correctly populates the underlying CDSs"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsInstance(gbk.genes[0], CDS)

    def test_populate_cds_mibig(self):
        """Tests whether parsing a mibig GBK correctly populates the underlying CDSs"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002609.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG, run)

        self.assertIsInstance(gbk.genes[0], CDS)

    def test_populate_hierarchical_objects(self):
        """Tests whether parsing a GBK correclty generates parent-child feature relations
        via checking for presence of the lowest level child - proto_core"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        proto_core = gbk.region.cand_clusters[1].proto_clusters[1].proto_core[1]

        self.assertIsInstance(proto_core, ProtoCore)

    def test_parse_gbk_has_dna_seq(self):
        """Tests whether parsing a GBK correclty has DNA sequence"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        dna_sequence = gbk.nt_seq

        self.assertIsInstance(dna_sequence, Seq)

    def test_filter_overlap_over_threshold(self):
        """Test whether add_cds_filter_overlap correctly throws out a cds"""

        gbk = GBK("", "", "")

        cds_a = CDS(0, 18)
        cds_a.aa_seq = "M" * 6
        cds_a.strand = 1

        cds_b = CDS(0, 9)
        cds_b.strand = 1
        cds_b.aa_seq = "M" * 3
        # nt_overlap_len_a_b = 9
        # aa_overlap = 9/3 = 3
        # 10% cds_b aa len = 0.1 * 3 = 0.3
        # aa_overlap > 10% shortest cds: 3 > 0.3

        gbk.add_cds_overlap_filter(cds_a, 0.1)

        gbk.add_cds_overlap_filter(cds_b, 0.1)

        expected_cds_list = [cds_a]

        self.assertEqual(expected_cds_list, gbk.genes)

    def test_filter_overlap_under_threshold(self):
        """Test whether filter_overlap correctly includes a new CDS if it overlaps with
        another CDS but under the cutoff threshold
        """

        gbk = GBK("", "", "")

        cds_a = CDS(0, 18)
        cds_a.aa_seq = "M" * 6
        cds_b = CDS(18, 36)
        cds_b.aa_seq = "M" * 4
        # nt_overlap_len_a_b = 1
        # aa_overlap = 1/3 = 0.33
        # 10% cds_b aa len = 0.1*4 = 0.4
        # aa_overlap < 10% shortest cds: 0.33 < 0.4

        gbk.add_cds_overlap_filter(cds_a, cds_overlap_cutoff=0.1)
        gbk.add_cds_overlap_filter(cds_b, cds_overlap_cutoff=0.1)

        expected_cds_list = [cds_a, cds_b]

        self.assertEqual(expected_cds_list, gbk.genes)

    def test_filter_overlap_over_threshold_diff_strands(self):
        """Test whether filter_overlap correclty preserves a CDS that is over the cds
        cutoff threshold, but is on a different strand
        """

        gbk = GBK("", "", "")

        cds_a = CDS(0, 18)
        cds_a.aa_seq = "M" * 6
        cds_a.strand = 1
        cds_b = CDS(0, 9)
        cds_b.aa_seq = "M" * 3
        cds_b.strand = -1
        # nt_overlap_len_a_b = 9
        # aa_overlap = 9/3 = 3
        # 10% cds_b aa len = 0.1 * 3 = 0.3
        # aa_overlap > 10% shortest cds: 3 > 0.3

        gbk.add_cds_overlap_filter(cds_a, cds_overlap_cutoff=0.1)
        gbk.add_cds_overlap_filter(cds_b, cds_overlap_cutoff=0.1)

        expected_cds_list = [cds_a, cds_b]

        self.assertEqual(expected_cds_list, gbk.genes)

    def test_save(self):
        """Tests whether a GBK object is correctly stored in the SQLite database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        gbk.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM gbk;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_save_all(self):
        """Tests whether this gbk and its children can all be saved to a database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        gbk.save_all()

        DB.commit()

        DB.save_to_disk(Path("tmp/db.db"))

        # 1 gbk, 11 bgc records
        expected_row_count = 12

        actual_row_count = 0

        # get gbk rows
        cursor_result = DB.execute_raw_query("SELECT * FROM gbk;")
        actual_row_count += len(cursor_result.fetchall())

        cursor_result = DB.execute_raw_query("SELECT * FROM bgc_record;")
        actual_row_count += len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_load_all(self):
        """Tests whether a set of GBKs can be recreated from a database"""
        populated_db_path = Path("test/test_data/database/valid_populated.db")
        DB.load_from_disk(populated_db_path)

        expected_gbk_count = 10

        actual_gbk_count = len(GBK.load_all())

        self.assertEqual(expected_gbk_count, actual_gbk_count)

    def test_load_many(self):
        """Tests whether a set of GBKs can be recreated from a database
        using load many
        """
        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY, run)

        gbk.save()

        loaded_gbks = GBK.load_many([gbk])

        self.assertEqual(gbk, loaded_gbks[0])

    def test_load_all_has_regions(self):
        """Tests whether the region objects were correctly loaded when loading GBKs"""
        populated_db_path = Path("test/test_data/database/valid_populated.db")
        DB.load_from_disk(populated_db_path)

        all_gbk = GBK.load_all()

        gbks_have_regions = [gbk.region is not None for gbk in all_gbk]

        self.assertTrue(all(gbks_have_regions))

    # TODO: test load candidate clusters, protoclusters, protocores

    def test_force_gbk_import(self):
        """Tests whether a GBK without any antiSMASH annotations can be imported
        when the --force-gbk option is enabled"""

        test_file = Path("test/test_data/invalid_gbk_folder/as4_no_cluster_feature.gbk")

        run = {
            "input_dir": Path("test/test_data/invalid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "force_gbk": True,
            "include_categories": set(),
            "exclude_categories": set(),
        }

        gbk = GBK.parse(test_file, SOURCE_TYPE.QUERY, run)

        actual_len = gbk.region.nt_stop
        expected_len = 31585

        self.assertEqual(expected_len, actual_len)
