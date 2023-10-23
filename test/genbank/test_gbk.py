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

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsInstance(gbk, GBK)

    def test_parse_as4gbk(self):
        """Tests whether an as4 GBK is instantiated correctly"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/CM001015.1.cluster001.gbk"
        )

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsInstance(gbk, GBK)

    def test_parse_as4_no_cluster_feature(self):
        """Tests whether an as4 gbk has no cluster feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as4_no_cluster_feature.gbk"
        )
        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY)

    def test_parse_metagenome_gbk(self):
        """Tests whether a metagenome GBK is instantiated correclty"""

        gbk_file_path = Path(
            "test/test_data/metagenome_valid_gbk_input/as5_metagenome_valid...region001.gbk"
        )

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsInstance(gbk, GBK)

    def test_parse_mibig_bac_region_only_gbk(self):
        """Tests whether a MIBiG bacterial gbk file, with only regions and CDSs, is instantiated correclty"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002476.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG)

        self.assertIsInstance(gbk, GBK)

    def test_parse_mibig_fun_region_only_gbk(self):
        """Tests whether a MIBiG fungal gbk file, with only regions and CDSs, is instantiated correclty"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002609.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG)

        self.assertIsInstance(gbk, GBK)

    def test_parse_gbk_multiple_regions(self):
        """Tests whether a GBK file has more than one region"""

        gbk_file_path = Path(
            "test/test_data/valid_gbk_multiple_regions_folder/valid_input_multiple_regions.gbk"
        )

        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY)

    def test_parse_as5_no_region(self):
        """Tests whether a GBK file has no region feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_no_region_feature.gbk"
        )

        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY)

    def test_parse_as5_no_cand_cluster(self):
        """Tests whether a GBK file with no cand_cluster features gives warning"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_no_cand_cluster_feature.gbk"
        )

        with self.assertLogs(level=logging.INFO) as cm:
            logging.info("nonsense")
            GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        # cm.output a list of strings of all the logs
        str = "does contain an antiSMASH cand_cluster feature"
        warning = any(str in log for log in cm.output)

        self.assertEqual(warning, True)

    def test_parse_as5_missing_protocluster(self):
        """Tests whether a GBK file has missing protocluster feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_missing_protocluster_feature.gbk"
        )

        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY)

    def test_parse_as5_missing_proto_core(self):
        """Tests whether a GBK file has missing proto_core feature"""
        gbk_file_path = Path(
            "test/test_data/invalid_gbk_folder/as5_missing_proto_core_feature.gbk"
        )

        self.assertRaises(InvalidGBKError, GBK.parse, gbk_file_path, SOURCE_TYPE.QUERY)

    def test_populate_regions(self):
        """Tests whether parsing a GBK correctly populates the underlying region"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsInstance(gbk.region, Region)

    def test_populate_cds(self):
        """Tests whether parsing a GBK correctly populates the underlying CDSs"""

        # GBK has one region
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        self.assertIsInstance(gbk.genes[0], CDS)

    def test_populate_cds_mibig(self):
        """Tests whether parsing a mibig GBK correctly populates the underlying CDSs"""

        gbk_file_path = Path("test/test_data/MIBiG_gbk/BGC0002609.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.MIBIG)

        self.assertIsInstance(gbk.genes[0], CDS)

    def test_populate_hierarchical_objects(self):
        """Tests whether parsing a GBK correclty generates parent-child feature relations
        via checking for presence of the lowest level child - proto_core"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        proto_core = gbk.region.cand_clusters[1].proto_clusters[1].proto_core[1]

        self.assertIsInstance(proto_core, ProtoCore)

    def test_parse_gbk_has_dna_seq(self):
        """Tests whether parsing a GBK correclty has DNA sequence"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        dna_sequence = gbk.nt_seq

        self.assertIsInstance(dna_sequence, Seq)

    def test_filter_overlap_over_threshold(self):
        """Test whether add_cds_filter_overlap correctly throws out a cds"""

        gbk = GBK("", "")

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

        gbk = GBK("", "")

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

        gbk = GBK("", "")

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

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

        gbk.save()

        cursor_result = DB.execute_raw_query("SELECT * FROM gbk;")

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_save_all(self):
        """Tests whether this gbk and its children can all be saved to a database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, SOURCE_TYPE.QUERY)

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

    def test_load_all_has_regions(self):
        """Tests whether the region objects were correctly loaded when loading GBKs"""
        populated_db_path = Path("test/test_data/database/valid_populated.db")
        DB.load_from_disk(populated_db_path)

        all_gbk = GBK.load_all()

        gbks_have_regions = [gbk.region is not None for gbk in all_gbk]

        self.assertTrue(all(gbks_have_regions))

    # TODO: test load candidate clusters, protoclusters, protocores
