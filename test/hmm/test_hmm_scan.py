"""Contains tests to test the hmm scanning functionality using pyhmmer"""

# from python
from pathlib import Path
import unittest
from unittest import TestCase
import platform

# from dependencies
from pyhmmer.easel import TextSequence

# from other modules
from big_scape.data import DB
from big_scape.hmm import HMMer, cds_to_input_task
from big_scape.genbank import GBK, CDS
from big_scape.hmm import HSP

import big_scape.enums as bs_enums


class TestHMMScan(TestCase):
    """Contains tests to check the hmmscan functionality"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_cds_to_input_task(self):
        """Tests the cds_to_input_task function, which generates input tasks for workers
        from cds objects"""
        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        cds = CDS(0, len(aa_seq) * 3)
        cds.aa_seq = aa_seq
        cds.orf_num = 1

        expected_result = (0, cds.aa_seq)
        actual_result = next(cds_to_input_task([cds]))

        self.assertEqual(expected_result, actual_result)

    def test_profile_hmmsearch(self):
        """Tests the profile_hmmsearch, which performs scanning of profiles for sequence
        domains"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        sequences = [TextSequence(name=str(0).encode(), sequence=aa_seq)]

        results = HMMer.profile_hmmsearch(sequences)

        expected_result = (0, "PF00457.19", 249.32315063476562, 60, 238)

        actual_result = results[0]

        self.assertEqual(expected_result, actual_result)

    def test_search(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        cds = CDS(0, len(aa_seq) * 3)
        cds.aa_seq = aa_seq

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        expected_result = HSP(cds, "PF00457", 249.32315063476562, 0, 0)

        HMMer.init(hmm_path)

        HMMer.hmmsearch_simple([cds])

        actual_result = cds.hsps[0]

        self.assertEqual(expected_result, actual_result)

    def test_search_two_domain_hits(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "VSMTLNGGGSYSTQWTNCGNFVAGKSVEGTKTFQVRYNGYFNPSGNGYGCLYGWTHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSGMQTYMSKNKCHCPNEVTIVSLSCHRTRTWWWVWEIPLNEDLVVGVVARSLYSMQQD"
            "GTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTMQQDGTQQDR"
            "IKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
        )

        # domain_hit_1
        # score: 113.76690673828125
        # included: True
        # env_start: 60
        # env_stop: 155

        # domain_hit_2
        # score: 40.309471130371094
        # included: True
        # env_start: 147
        # env_stop: 204

        # domain_hit_3
        # score: 37.029972076416016
        # included: True
        # env_start: 203
        # env_stop: 238

        # domain_hit_4
        # score: -1.5011663436889648
        # included: False
        # env_start: 0
        # env_stop: 95

        # domain_hit_5
        # score: 0.8793293833732605
        # included: False
        # env_start: 0
        # env_stop: 95

        # domain_hit_1 and domain_hit_2 overlap, so we only expect domain_hit_1 to be included,
        # as well as the non-overlapping domain_hit_3

        cds = CDS(0, len(aa_seq) * 3)
        cds.aa_seq = aa_seq

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        HMMer.hmmsearch_simple([cds])

        expected_result = [113.76690673828125, 37.029972076416016]
        actual_result = [domain.score for domain in cds.hsps]

        self.assertEqual(expected_result, actual_result)

    @unittest.skipIf(platform.system() == "Darwin", "Unsupported OS")
    def test_scan(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        cds = CDS(0, len(aa_seq) * 3)
        cds.aa_seq = aa_seq

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        expected_result = HSP(cds, "PF00457", 249.32315063476562, 0, 0)

        HMMer.hmmsearch_multiprocess([cds], batch_size=1)

        actual_result = cds.hsps[0]

        self.assertEqual(expected_result, actual_result)

    @unittest.skipIf(platform.system() == "Darwin", "Unsupported OS")
    def test_search_two_domain_hits_multiprocess(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "VSMTLNGGGSYSTQWTNCGNFVAGKSVEGTKTFQVRYNGYFNPSGNGYGCLYGWTHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSGMQTYMSKNKCHCPNEVTIVSLSCHRTRTWWWVWEIPLNEDLVVGVVARSLYSMQQD"
            "GTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTMQQDGTQQDR"
            "IKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
        )

        cds = CDS(0, len(aa_seq) * 3)
        cds.aa_seq = aa_seq

        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        HMMer.hmmsearch_multiprocess([cds], batch_size=1)

        expected_result = [113.76690673828125, 37.029972076416016]
        actual_result = [domain.score for domain in cds.hsps]

        self.assertEqual(expected_result, actual_result)

    def test_save_hsp(self):
        """Tests whether an HSP can be correctly saved to the database"""
        DB.create_in_mem()

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        gbk = GBK("", "", bs_enums.SOURCE_TYPE.QUERY)
        gbk.metadata = {"organism": "test", "taxonomy": "test", "description": "test"}
        gbk.save()

        cds = CDS(0, len(aa_seq) * 3)
        cds.strand = 1
        cds.orf_num = 1
        cds.aa_seq = aa_seq
        cds.parent_gbk = gbk
        cds.save()

        hsp = HSP(cds, "PF00457", 249.32315063476562, 0, 0)
        hsp.save()

        expected_row_count = 1

        actual_row_count = 0

        cursor_result = DB.execute_raw_query("SELECT * FROM hsp;")
        actual_row_count += len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
