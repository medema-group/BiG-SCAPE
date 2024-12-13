"""Contains tests for the BGCRecord class, a generalized class from which other genbank
clasess inherit
"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS, Region
from big_scape.hmm import HSP
from big_scape.utility import domain_includelist_filter
from big_scape import enums as bs_enums
from big_scape.file_input.load_files import get_all_bgc_records


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
        actual_cds_count = len(list(record.get_cds()))

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
        actual_cds_count = len(list(record.get_cds(True)))

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

    def test_domainincludelist_filter_all(self):
        """Tests whether the domain_includelist_filter all
        correctly filters out records that do not contain all domains"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "record_type": bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER,
            "domain_includelist_all": ["PF00001", "PF00002"],
            "domain_includelist_any": None,
        }

        gbk_1 = GBK("", "", "")
        gbk_1.region = Region(gbk_1, 0, 0, 100, False, "")
        cds_1 = CDS(10, 90)
        cds_1.strand = 1
        gbk_1.genes.append(cds_1)
        domains_1 = ["PF00001", "PF00002", "PF00003", "PF00004", "PF00005"]
        for domain in domains_1:
            cds_1.hsps.append(HSP(cds_1, domain, 100, 0, 30))

        gbk_2 = GBK("", "", "")
        gbk_2.region = Region(gbk_2, 0, 0, 100, False, "")
        cds_2 = CDS(10, 90)
        cds_2.strand = 1
        gbk_2.genes.append(cds_2)
        domains_2 = ["PF00001", "PF00003", "PF00004", "PF00005"]
        for domain in domains_2:
            cds_2.hsps.append(HSP(cds_2, domain, 100, 0, 30))

        gbks = [gbk_1, gbk_2]

        all_bgc_records = get_all_bgc_records(run, gbks)

        domainlist_bgc_records = domain_includelist_filter(run, all_bgc_records)

        expected_records = get_all_bgc_records(run, [gbk_1])

        self.assertEqual(expected_records, domainlist_bgc_records)

    def test_domainincludelist_filter_any(self):
        """Tests whether the domain_includelist_filter any
        correctly filters out records that do not contain all domains"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "record_type": bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER,
            "domain_includelist_all": None,
            "domain_includelist_any": ["PF00001", "PF00002"],
        }

        gbk_1 = GBK("", "", "")
        gbk_1.region = Region(gbk_1, 0, 0, 100, False, "")
        cds_1 = CDS(10, 90)
        cds_1.strand = 1
        gbk_1.genes.append(cds_1)
        domains_1 = ["PF00002", "PF00003", "PF00004", "PF00005"]
        for domain in domains_1:
            cds_1.hsps.append(HSP(cds_1, domain, 100, 0, 30))

        gbk_2 = GBK("", "", "")
        gbk_2.region = Region(gbk_2, 0, 0, 100, False, "")
        cds_2 = CDS(10, 90)
        cds_2.strand = 1
        gbk_2.genes.append(cds_2)
        domains_2 = ["PF00001", "PF00003", "PF00004", "PF00005"]
        for domain in domains_2:
            cds_2.hsps.append(HSP(cds_2, domain, 100, 0, 30))

        gbks = [gbk_1, gbk_2]

        all_bgc_records = get_all_bgc_records(run, gbks)

        domainlist_bgc_records = domain_includelist_filter(run, all_bgc_records)

        self.assertEqual(all_bgc_records, domainlist_bgc_records)

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
        actual_cds_count = len(list(region.get_cds_with_domains()))

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
