"""Contains tests for the CDS class and functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from Bio.SeqFeature import (
    SeqFeature,
    FeatureLocation,
    Seq,
    BeforePosition,
    AfterPosition,
)

# from other modules
from src.hmm import HSP
from src.genbank import (
    GBK,
    CDS,
    check_translation,
    get_translation,
    translate,
    trim_fuzzy,
)
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

    def test_parse_aa_seq(self):
        """Tests whether an aa sequence is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["gene_kind"] = "biosynthetic"

        expected_aa_seq = "MDRAAGMSVRLK"

        feature.qualifiers["translation"] = [expected_aa_seq]

        parent_gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region.gbk"
        )
        parent_gbk = GBK.parse(parent_gbk_file_path, "query")

        cds = CDS.parse(feature, parent_gbk)

        self.assertEqual(expected_aa_seq, cds.aa_seq)

    def test_parse_gene_kind(self):
        """Tests whether a gene_kind is correclty parsed from a feature"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        feature.qualifiers["translation"] = ["MDRAAGMSVRLK"]

        expected_gene_kind = "biosynthetic"

        feature.qualifiers["gene_kind"] = [expected_gene_kind]

        parent_gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region.gbk"
        )
        parent_gbk = GBK.parse(parent_gbk_file_path, "query")

        cds = CDS.parse(feature, parent_gbk)

        self.assertEqual(expected_gene_kind, cds.gene_kind)

    def test_parse_wrong_type(self):
        """Tests whether create_region correctly throws an error when given a feature of
        a wrong type
        """

        feature = SeqFeature(type="region")

        parent_gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region.gbk"
        )
        parent_gbk = GBK.parse(parent_gbk_file_path, "query")

        self.assertRaises(InvalidGBKError, CDS.parse, feature, parent_gbk)

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

    def test_load_all(self):
        """Tests whether a set of CDS objects can be recreated from a database"""
        populated_db_path = Path("test/test_data/database/valid_populated.db")
        DB.load_from_disk(populated_db_path)

        expected_cds_count = 355

        all_gbk = GBK.load_all()

        gbk_cds_count = [len(gbk.genes) for gbk in all_gbk]

        actual_gbk_count = sum(gbk_cds_count)

        self.assertEqual(expected_cds_count, actual_gbk_count)

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

    def test_len_overlap(self):
        """Tests whether the len_overlap function returns a correct overlap len"""

        cds_a = CDS(0, 50)
        cds_b = CDS(10, 80)

        expected_len = 40

        actual_len = CDS.len_nt_overlap(cds_a, cds_b)

        self.assertEqual(expected_len, actual_len)

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

    def test_translate_output_type(self):
        """Tests whether translate correclty outputs a Seq object"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("ATGCAGCAGGACGGCACACAGCAGGACCGGATCAAGCAGAGTCCCGCCCCTCTCTGA")

        transl_nt_seq = translate(feature, nt_seq)

        self.assertIsInstance(transl_nt_seq, Seq)

    def test_translate_correct_output(self):
        """Tests whether translate translates correclty"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("ATGCAGCAGGACGGCACACAGCAGGACCGGATCAAGCAGAGTCCCGCCCCTCTCTGA")

        expected_transl_nt_seq = Seq("MQQDGTQQDRIKQSPAPL")
        transl_nt_seq = translate(feature, nt_seq)

        self.assertEqual(expected_transl_nt_seq, transl_nt_seq)

    def test_check_translation_true(self):
        """Tests whether check_translation returns a correct assessment"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("CAGCAGGACGGCACACAGCAGGACCGGATCAAGCAGAGTCCCGCCCCTCTCTGA")
        aa_seq = Seq("QQDGTQQDRIKQSPAPL")

        expected_assessment = True
        assessment = check_translation(aa_seq, nt_seq, feature)

        self.assertEqual(expected_assessment, assessment)

    def test_trim_fuzzy_with_fuzzy_start(self):
        """Test of trim_fuzzy trims a fuzzy start correclty"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("ATGCAGCAGGACGGCACAC")
        fuzzy_start = True
        fuzzy_end = False
        remainder = 2

        expected_nt_seq = Seq("GCAGCAGGACGGCACAC")

        trimmed_nt_seq = trim_fuzzy(feature, nt_seq, fuzzy_start, fuzzy_end, remainder)

        self.assertEqual(expected_nt_seq, trimmed_nt_seq)

    def test_trim_fuzzy_with_fuzzy_end(self):
        """Test of trim_fuzzy trims a fuzzy end correclty"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("ATGCAGCAGGACGGCACAC")
        fuzzy_start = False
        fuzzy_end = True
        remainder = 2

        expected_nt_seq = Seq("ATGCAGCAGGACGGCAC")

        trimmed_nt_seq = trim_fuzzy(feature, nt_seq, fuzzy_start, fuzzy_end, remainder)

        self.assertEqual(expected_nt_seq, trimmed_nt_seq)

    def test_trim_fuzzy_error_both_fuzzies(self):
        """Test if trim_fuzzy throws an error when passing the same value for
        both fuzzies"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")
        nt_seq = Seq("ATGCAGCAGGACGGCACAC")
        fuzzy_start = True
        fuzzy_end = True
        remainder = 2

        self.assertRaises(
            ValueError, trim_fuzzy, feature, nt_seq, fuzzy_start, fuzzy_end, remainder
        )

    def test_get_translation_same_fuzzies_true(self):
        """Test if get_translation returns None when given the same value (true) for both fuzzies"""

        feature = SeqFeature(
            FeatureLocation(BeforePosition(5), AfterPosition(10), strand=1), type="CDS"
        )

        nt_seq = Seq("ATGCAGCAGGACGGCACAC")

        expected_translation = None

        translation = get_translation(feature, nt_seq)

        self.assertEqual(expected_translation, translation)

    def test_get_translation_same_fuzzies_false(self):
        """Test if get_translation returns None when given the same value (false) for both fuzzies"""

        feature = SeqFeature(FeatureLocation(5, 10, strand=1), type="CDS")

        nt_seq = Seq("ATGCAGCAGGACGGCACAC")

        expected_translation = None

        translation = get_translation(feature, nt_seq)

        self.assertEqual(expected_translation, translation)

    def test_get_translation_return_seq_object(self):
        """Test if get_translation returns None when given the same value (false) for both fuzzies"""

        feature = SeqFeature(
            FeatureLocation(BeforePosition(5), 10, strand=1), type="CDS"
        )

        nt_seq = Seq("ATGCAGCAGGACGGCACAC")

        translation = get_translation(feature, nt_seq)

        self.assertIsInstance(translation, Seq)

    def test_parse_translate_aa_seq_cds_1(self):
        """Tests whether an aa sequence is correclty translated from a feature
        with no fuzzies: complement(770..1447), and no translation available"""

        parent_gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region_cds_no_trans.gbk"
        )
        parent_gbk = GBK.parse(parent_gbk_file_path, "query")
        cds_1 = parent_gbk.genes[1]

        expected_translation = (
            "MRVLIVEDEPYLAEAIRDGLRLEAIAADTAGNGDTALELLSLNTY"
            "DIAVLDRDIPGPSGDEIAKRIVASGSGLPILMLTAADRLDDKITGFELGADDYLTKPFE"
            "LRELVLRLRALDRRRAHNRPPVLEIAGLRLNPFRREVYRDDRYIALTRKQFAVLEVLVS"
            "ADGGVVSAEELLERAWDKNADPFTNAVRITVSALRKRLGEPWIITTVAGVGYRIGAAPG"
            "AGR"
        )

        translation = cds_1.aa_seq

        self.assertEqual(expected_translation, translation)

    def test_parse_translate_aa_seq_cds_0(self):
        """Tests whether an aa sequence is correclty translated from a feature
        with location: complement(<1..759) and no translation available"""

        parent_gbk_file_path = Path(
            "test/test_data/valid_gbk_folder/valid_input_region_cds_no_trans.gbk"
        )
        parent_gbk = GBK.parse(parent_gbk_file_path, "query")
        cds_0 = parent_gbk.genes[0]

        expected_translation = (
            "MDRAAGMSVRLKLTLSYACFLVLAGVLLLASVWLFLLRDVPDVLA"
            "KPPPGGVLERSVLVRNFLPAAGSVLFFLLLFGLLGGWILAGRMLAPLTRITDAARMAAN"
            "GSLSHRIRLEGTEDEFRELADAFDAMLARLEAHVAAQRRFAANASHELRTPLAITQALL"
            "EVARNDPAKDPLLVFDRLHAVNARAIDLTEALLVLSRADQRAFTREPVDLSLLVEEAIE"
            "TLLPIAEKRRVVIIASGHISRVVGSATLLLQ"
        )

        translation = cds_0.aa_seq

        self.assertEqual(expected_translation, translation)
