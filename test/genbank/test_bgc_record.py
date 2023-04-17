"""Contains tests for the BGCRecord class, a generalized class from which other genbank
clasess inherit
"""

# from python
from unittest import TestCase

# from other modules
from src.genbank import GBK, BGCRecord, CDS


class TestBGCRecord(TestCase):
    """Contains tests for the BGC record class"""

    def test_get_cds(self):
        """Tests whether the get_cds method on this BGCRecord class correctly retrieves
        the CDS objects in the parent GBK that lie within the range of the record's
        nt_start and nt_stop coordinates
        """

        gbk = GBK("", "")
        record = BGCRecord()
        record.parent_gbk = gbk
        record.nt_start = 10
        record.nt_stop = 90

        gbk.genes = [
            CDS(0, 10),
            CDS(10, 20),
            CDS(20, 30),
            CDS(30, 40),
            CDS(40, 50),
            CDS(50, 60),
            CDS(60, 70),
            CDS(70, 80),
            CDS(80, 90),
            CDS(90, 100),
        ]  # 10 total

        # coordinates are exclusive. first and last in above list should not be included
        expected_cds_count = 8

        actual_cds_count = len(record.get_cds())

        self.assertEqual(expected_cds_count, actual_cds_count)
