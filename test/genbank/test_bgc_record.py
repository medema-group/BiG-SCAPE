"""Contains tests for the BGCRecord class, a generalized class from which other genbank
clasess inherit
"""

# from python
from unittest import TestCase

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS, Region
from big_scape.hmm import HSP


class TestBGCRecord(TestCase):
    """Contains tests for the BGC record class"""

    def test_get_cds(self):
        """Tests whether the get_cds method on this BGCRecord class correctly retrieves
        the CDS objects in the parent GBK that lie within the range of the record's
        nt_start and nt_stop coordinates
        """

        gbk = GBK("", "", "")
        record = BGCRecord(gbk, 0, 10, 90, False, "")

        # 10 total
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
        ]

        # coordinates are exclusive. first and last in above list should not be included
        expected_cds_count = 8

        # TODO: use __len__
        actual_cds_count = len(record.get_cds())

        self.assertEqual(expected_cds_count, actual_cds_count)

    def test_get_cds_all(self):
        """Tests whether the get_cds method on this BGCRecord class correctly retrieves
        all CDS objects in the parent GBK when get_cds is called with the return_all
        parameter set to True
        """

        gbk = GBK("", "", "")
        record = BGCRecord(gbk, 0, 10, 90, False, "")

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
        ]

        # coordinates are exclusive. first and last in above list should not be included
        expected_cds_count = 10

        # TODO: use __len__
        actual_cds_count = len(record.get_cds(True))

        self.assertEqual(expected_cds_count, actual_cds_count)

    def test_get_hsps(self):
        """Tests whether the get_hsps method on this BGCRecord class correctly retreives
        the HPS objects that are in range of all CDSes that belong to this record in a
        GBK
        """
        domains = ["PF00001", "PF00002", "PF00003", "PF00004", "PF00005"]

        gbk = GBK("", "", "")
        gbk.region = Region(gbk, 0, 0, 100, False, "")
        cds = CDS(10, 90)
        cds.strand = 1
        gbk.genes.append(cds)

        for domain in domains:
            cds.hsps.append(HSP(cds, domain, 100, 0, 30))

        expected_domains = domains

        hsps = gbk.region.get_hsps()

        actual_domains = [hsp.domain for hsp in hsps]

        self.assertEqual(expected_domains, actual_domains)

    def test_get_cds_with_domains(self):
        """Tests whether the test_get_cds_with_domains method correctly retrieves a
        subset of CDS containing only domains
        """
        cds_list = [
            CDS(0, 100),
            CDS(0, 100),
            CDS(0, 100),
            CDS(0, 100),
            CDS(0, 100),
            CDS(0, 100),
            CDS(0, 100),
        ]

        cds_list[0].hsps = [HSP(cds_list[0], "test", 1.0, 0, 100)]
        cds_list[2].hsps = [HSP(cds_list[0], "test", 1.0, 0, 100)]
        cds_list[4].hsps = [HSP(cds_list[0], "test", 1.0, 0, 100)]
        cds_list[6].hsps = [HSP(cds_list[0], "test", 1.0, 0, 100)]

        gbk = GBK("", "", "test")
        gbk.genes = cds_list
        region = BGCRecord(gbk, 0, 0, 100, False, "")

        expected_cds_count = 4
        actual_cds_count = len(region.get_cds_with_domains())

        self.assertEqual(expected_cds_count, actual_cds_count)

    def test_get_cds_start_stop_region(self):
        """Tests whether get_cds_start_stop correctly returns full region"""
        gbk = GBK("", "", "test")
        cds = [
            CDS(0, 20),
            CDS(20, 30),
            CDS(30, 40),
            CDS(40, 50),
            CDS(50, 60),
            CDS(80, 100),
        ]
        gbk.genes = cds
        region = BGCRecord(gbk, 0, 0, 100, False, "")

        expected_slice = (1, 6)
        actual_slice = region.get_cds_start_stop()

        self.assertEqual(expected_slice, actual_slice)

    def test_get_cds_start_stop_protocluster(self):
        """Tests whether get_cds_start_stop correctly returns full region"""
        gbk = GBK("", "", "test")
        cds = [
            CDS(0, 20),
            CDS(20, 30),
            CDS(30, 40),
            CDS(40, 50),
            CDS(50, 60),
            CDS(80, 100),
        ]
        gbk.genes = cds
        protocluster = BGCRecord(gbk, 0, 30, 60, False, "")

        expected_slice = (3, 5)
        actual_slice = protocluster.get_cds_start_stop()

        self.assertEqual(expected_slice, actual_slice)

    def test_get_cds_start_stop_partial_region_end(self):
        """Tests partial region bounds"""
        gbk = GBK("", "", "test")
        cds = [
            CDS(0, 20),
            CDS(20, 30),
            CDS(30, 40),
            CDS(40, 50),
            CDS(50, 60),
            CDS(80, 100),
        ]
        gbk.genes = cds
        region = BGCRecord(gbk, 0, 0, 55, False, "")

        expected_slice = (1, 4)
        actual_slice = region.get_cds_start_stop()

        self.assertEqual(expected_slice, actual_slice)

    def test_get_cds_start_stop_partial_region_start(self):
        """Tests partial region bounds"""
        gbk = GBK("", "", "test")
        cds = [
            CDS(0, 20),
            CDS(20, 30),
            CDS(30, 40),
            CDS(40, 50),
            CDS(50, 60),
            CDS(80, 100),
        ]
        gbk.genes = cds
        region = BGCRecord(gbk, 0, 5, 100, False, "")

        expected_slice = (2, 6)
        actual_slice = region.get_cds_start_stop()

        self.assertEqual(expected_slice, actual_slice)
