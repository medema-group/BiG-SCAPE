"""Contains tests for the CDS class and functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
from src.genbank import GBK, CDS
from src.hmm import HSP
from src.errors import InvalidGBKError
from src.data import DB


class TestCDS(TestCase):
    """Test class for base GBK parsing tests"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_create_cds(self):
        """Tests whether a cds is instatiated correclty"""

        expected_start = 1
        expected_end = 5

        cds = CDS(expected_start, expected_end)

        self.assertIsInstance(cds, CDS)

    def test_parse_aa_seq(self):
        """Tests whether an aa sequence is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["gene_kind"] = "biosynthetic"

        expected_aa_seq = "MDRAAGMSVRLK"

        feature.qualifiers["translation"] = [expected_aa_seq]

        cds = CDS.parse(feature)

        self.assertEqual(expected_aa_seq, cds.aa_seq)

    def test_parse_no_aa_seq(self):
        """Tests whehter parse correctly throws an error when given a feature
        lacking a translation qualifier
        """

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["gene_kind"] = ["biosynthetic"]

        self.assertRaises(InvalidGBKError, CDS.parse, feature)

    def test_parse_gene_kind(self):
        """Tests whether a gene_kind is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
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

    def test_save_cds(self):
        """Tests whether a CDS can be correctly saved to the database"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")

        gbk = GBK.parse(gbk_file_path, "query")

        gbk.save_all()

        DB.commit()

        DB.save_to_disk(Path("tmp/db.db"))

        # 1 gbk, 37 bgc records
        expected_row_count = 37

        actual_row_count = 0

        cursor_result = DB.execute_raw_query("SELECT * FROM cds;")
        actual_row_count += len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

    def test_cds_has_overlap_false(self):
        """Tests the has_overlap function where a is left of b"""

        cds_a = CDS(0, 50)
        cds_b = CDS(100, 150)

        expected_result = False

        actual_result = CDS.has_overlap(cds_a, cds_b)

        self.assertEqual(expected_result, actual_result)

    def test_cds_has_overlap_true(self):
        """Tests the has_overlap function where a is left of b"""

        cds_a = CDS(0, 50)
        cds_b = CDS(10, 80)

        expected_result = True

        actual_result = CDS.has_overlap(cds_a, cds_b)

        self.assertEqual(expected_result, actual_result)

    def test_cds_has_overlap_diff_strands(self):
        """Tests the has_overlap function where a and b are in different strands"""

        cds_a = CDS(0, 50)
        cds_a.strand = 1
        cds_b = CDS(10, 80)
        cds_b.strand = -1

        expected_result = False

        actual_result = CDS.has_overlap(cds_a, cds_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap(self):
        """Tests whether the len_overlap function returns a correct overlap len"""

        cds_a = CDS(0, 50)
        cds_b = CDS(10, 80)

        expected_len = 40

        actual_len = CDS.len_nt_overlap(cds_a, cds_b)

        self.assertEqual(expected_len, actual_len)

    def test_len_overlap_diff_strands(self):
        """Tests whether the len_overlap function returns len = 0 if cds in different
        strands"""

        cds_a = CDS(0, 50)
        cds_a.strand = 1
        cds_b = CDS(10, 80)
        cds_b.strand = -1

        expected_result = 0

        actual_result = CDS.len_nt_overlap(cds_a, cds_b)

        self.assertEqual(expected_result, actual_result)

    def test_filter_overlap_over_threshold(self):
        """Test whether filter_overlap correclty throws out a cds"""

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

        expected_cds_list = [cds_a]

        cds_list = [cds_a, cds_b]

        CDS.filter_overlap(cds_list, 0.1)

        self.assertEqual(expected_cds_list, cds_list)

    def test_filter_overlap_under_threshold(self):
        """Test whether filter_overlap correclty throws out a cds"""

        cds_a = CDS(0, 18)
        cds_a.aa_seq = "M" * 6
        cds_b = CDS(18, 36)
        cds_b.aa_seq = "M" * 4
        # nt_overlap_len_a_b = 1
        # aa_overlap = 1/3 = 0.33
        # 10% cds_b aa len = 0.1*4 = 0.4
        # aa_overlap < 10% shortest cds: 0.33 < 0.4

        expected_cds_list = [cds_a, cds_b]

        cds_list = [cds_a, cds_b]

        CDS.filter_overlap(cds_list, 0.1)

        self.assertEqual(expected_cds_list, cds_list)

    def test_filter_overlap_over_threshold_diff_strands(self):
        """Test whether filter_overlap correclty throws out a cds"""

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

        expected_cds_list = [cds_a, cds_b]

        cds_list = [cds_a, cds_b]

        CDS.filter_overlap(cds_list, 0.1)

        self.assertEqual(expected_cds_list, cds_list)

    def test_add_hsp_overlap_filter(self):
        """Tests whether a new hsp can be added on a new CDS using the add with overlap
        filter function
        """
        cds = CDS(0, 100)
        cds.aa_seq = "M" * 100
        cds.strand = 1

        new_hsp = HSP(cds, "test_domain", 100.0, 0, 100)

        cds.add_hsp_overlap_filter(new_hsp)

        expected_result = [new_hsp]

        actual_result = cds.hsps

        self.assertEqual(expected_result, actual_result)

    def test_add_hsp_overlap_filter_keep_old(self):
        """Tests whether a newly added HSP is discarded because it overlaps with a
        better scoring existing HSP
        """
        cds = CDS(0, 100)
        cds.aa_seq = "M" * 100
        cds.strand = 1

        old_hsp = HSP(cds, "test_domain_1", 100.0, 0, 100)

        cds.add_hsp_overlap_filter(old_hsp)

        new_hsp = HSP(cds, "test_domain_2", 50.0, 0, 100)

        cds.add_hsp_overlap_filter(new_hsp)

        expected_result = [old_hsp]

        actual_result = cds.hsps

        self.assertEqual(expected_result, actual_result)

    def test_add_hsp_overlap_filter_keep_new(self):
        """Tests whether an old HSP is replaced because it overlaps with a better
        scoring new HSP
        """
        cds = CDS(0, 100)
        cds.aa_seq = "M" * 100
        cds.strand = 1

        old_hsp = HSP(cds, "test_domain_1", 50.0, 0, 100)

        cds.add_hsp_overlap_filter(old_hsp)

        new_hsp = HSP(cds, "test_domain_2", 100.0, 0, 100)

        cds.add_hsp_overlap_filter(new_hsp)

        expected_result = [new_hsp]

        actual_result = cds.hsps

        self.assertEqual(expected_result, actual_result)

    def test_add_hsp_overlap_filter_keep_new_same_bitscore(self):
        """Tests whether an old HSP is replaced because it has a higher env_start than
        a newly added HSP with the same bitscore
        """
        cds = CDS(0, 100)
        cds.aa_seq = "M" * 100
        cds.strand = 1

        old_hsp = HSP(cds, "test_domain_1", 100.0, 20, 100)

        cds.add_hsp_overlap_filter(old_hsp)

        new_hsp = HSP(cds, "test_domain_2", 100.0, 0, 100)

        cds.add_hsp_overlap_filter(new_hsp)

        expected_result = [new_hsp]

        actual_result = cds.hsps

        self.assertEqual(expected_result, actual_result)
