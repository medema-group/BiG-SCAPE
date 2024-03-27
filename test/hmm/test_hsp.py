"""Contains tests for the HSP class"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.genbank import GBK, CDS
from big_scape.data import DB
from big_scape.hmm import HSP


class TestHSP(TestCase):
    """Contains tests for HSP class functions"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_load_all(self):
        DB.load_from_disk(Path("test/test_data/database/valid_populated.db"))

        all_gbks = GBK.load_all()

        expected_hsp_count = 501
        actual_hsp_count = 0

        for gbk in all_gbks:
            HSP.load_all(gbk.genes)
            actual_hsp_count += sum([len(cds.hsps) for cds in gbk.genes])

        self.assertEqual(expected_hsp_count, actual_hsp_count)

    def test_equality(self):
        """Tests whether the equality operator works as expected for HSP objects, i.e.
        two HSP objects are equal if they have the same sanitized domain"""

        hsp1 = HSP(1, "PF01234.1", 1.0, 0, 10)
        hsp2 = HSP(1, "PF01234.2", 1.0, 0, 100)

        self.assertEqual(hsp1, hsp2)

    def test_inequality(self):
        """Tests whether the inequality operator works as expected for HSP objects, i.e.
        two HSP objects are not equal if they have different sanitized domains"""

        hsp1 = HSP(1, "PF01232.1", 1.0, 0, 100)
        hsp2 = HSP(1, "PF01234.2", 1.0, 0, 100)

        self.assertNotEqual(hsp1, hsp2)

    def test_equality_error_not_hsp(self):
        """Tests whether the equality operator works as expected for HSP objects, i.e.
        two HSP objects are not equal if the other object is not an HSP object"""

        hsp = HSP(1, "PF01234.1", 1.0, 0, 100)

        self.assertRaises(NotImplementedError, lambda: hsp == 1)

    def test_hsp_greater_than_env_based(self):
        """Tests whether the greater than operator works as expected for HSP objects, i.e.
        an HSP object is greater than another if:
        - self.cds.nt_start > other.cds.nt_start
        - self.env_start > other.env_start
        - !! self.env_end > other.env_end !!
        - self.score > other.score
        """

        cds_self = CDS(120, 100)
        cds_other = CDS(100, 150)

        hsp_self = HSP(cds_self, "PF01234.1", 1.0, 10, 100)
        hsp_other = HSP(cds_other, "PF01234.2", 2.0, 0, 50)

        val = hsp_self > hsp_other

        self.assertTrue(val)

    def test_hsp_greater_than_cds_start_lesser(self):
        """Tests whether the greater than operator works as expected for HSP objects, i.e.
        an HSP object is greater than another if:
        - !! self.cds.nt_start < other.cds.nt_start !!
        - self.env_start > other.env_start
        - self.env_end > other.env_end
        - self.score > other.score
        """

        cds_self = CDS(100, 100)
        cds_other = CDS(120, 150)

        hsp_self = HSP(cds_self, "PF01234.1", 1.0, 10, 100)
        hsp_other = HSP(cds_other, "PF01234.2", 2.0, 0, 50)

        val = hsp_self > hsp_other

        self.assertFalse(val)

    def test_hsp_greater_than_score_based(self):
        """Tests whether the greater than operator works as expected for HSP objects, i.e.
        an HSP object is greater than another if:
        - self.cds.nt_start > other.cds.nt_start
        - self.env_start = other.env_start
        - self.env_start = other.env_start
        - self.score > other.score
        """

        cds_self = CDS(120, 100)
        cds_other = CDS(100, 150)

        hsp_self = HSP(cds_self, "PF01234.1", 1.0, 10, 100)
        hsp_other = HSP(cds_other, "PF01234.2", 2.0, 10, 50)

        val = hsp_self > hsp_other

        self.assertTrue(val)

    def test_hsp_lesser_than_score_based(self):
        """Tests whether the greater than operator works as expected for HSP objects, i.e.
        an HSP object is greater than another if:
        - self.cds.nt_start > other.cds.nt_start
        - self.env_start = other.env_start
        - self.env_start = other.env_start
        - self.score > other.score
        """

        cds_self = CDS(120, 100)
        cds_other = CDS(100, 150)

        hsp_self = HSP(cds_self, "PF01234.1", 1.0, 10, 100)
        hsp_other = HSP(cds_other, "PF01234.2", 2.0, 10, 50)

        val = hsp_self > hsp_other

        self.assertTrue(val)

    def test_sorting_hsps(self):
        """Tests whether a list of HSP objects can be sorted correctly"""

        cds_self = CDS(120, 100)
        cds_other = CDS(100, 150)

        hsp_self = HSP(cds_self, "PF01234.1", 1.0, 10, 100)
        hsp_other = HSP(cds_other, "PF01234.2", 2.0, 10, 50)

        hsps = [hsp_other, hsp_self]
        hsps.sort()

        self.assertEqual(hsps, [hsp_self, hsp_other])

    def test_sanitize_accession(self):
        """tests the sanitize_accession method which strips the version from the pfam accession"""

        accession = "PF01234.1"

        accession_sanitized = HSP.sanitize_accession(accession)

        self.assertEqual(accession_sanitized, "PF01234")

    def test_sanitize_accession_no_version(self):
        """tests the sanitize_accession method which strips the version from the pfam accession"""

        accession = "PF01234"

        accession_sanitized = HSP.sanitize_accession(accession)

        self.assertEqual(accession_sanitized, "PF01234")

    def test_sanitize_accession_too_long(self):
        """Tests the sanitize_accession method which raises an error if the accession is too long"""

        accession = "PF0123453453"

        self.assertRaises(ValueError, lambda: HSP.sanitize_accession(accession))
