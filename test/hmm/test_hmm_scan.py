"""Contains tests to test the hmm scanning functionality using pyhmmer"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from pyhmmer.easel import TextSequence

# from other modules
from src.hmm import HMMer
from src.hmm.hmm import cds_to_input_task, has_overlap, len_overlap
from src.genbank import CDS
from src.hmm import HSP


class TestHMMScan(TestCase):
    """Contains tests to check the hmmscan functionality"""

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

        expected_result = (0, "PF00457.19", 249.32315063476562)

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

        expected_result = HSP(cds, "PF00457.19", 249.32315063476562)

        HMMer.init(hmm_path)

        actual_result = next(HMMer.search([cds]))

        self.assertEqual(expected_result, actual_result)

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

        expected_result = HSP(cds, "PF00457.19", 249.32315063476562)

        actual_result = next(HMMer.scan([cds]))

        self.assertEqual(expected_result, actual_result)

    def test_no_overlap_left(self):
        """Tests the no_overlap function where a is left of b"""
        cds_a = CDS(0, 50)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = False

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_no_overlap_right(self):
        """Tests the no_overlap function where a is right of b"""
        cds_a = CDS(100, 150)
        cds_b = CDS(0, 50)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = False

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_full(self):
        """Tests the no_overlap function for two overlapping regions

        both regions overlap perfectly"""
        cds_a = CDS(100, 150)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_partial_left(self):
        """Tests the no_overlap function for two overlapping regions

        a overlaps b on b left side"""
        cds_a = CDS(0, 120)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_has_overlap_partial_right(self):
        """Tests the no_overlap function for two overlapping regions

        a overlaps b on b right side"""
        cds_a = CDS(130, 250)
        cds_b = CDS(100, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = True

        actual_result = has_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_full(self):
        """Tests the len_overlap function for two fully overlapping regions"""
        cds_a = CDS(0, 100)
        cds_b = CDS(0, 100)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 100

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_left(self):
        """Tests the len_overlap function for two overlapping regions with 50 bp overlap

        a overlaps b on b left side"""
        cds_a = CDS(0, 100)
        cds_b = CDS(50, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 50

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_overlap_right(self):
        """Tests the len_overlap function for two overlapping regions with 50 bp overlap

        a overlaps b on b right side"""
        cds_a = CDS(100, 200)
        cds_b = CDS(50, 150)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 50

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_no_overlap_left(self):
        """Tests the len_overlap function for two non-overlapping regions

        a < b"""
        cds_a = CDS(0, 100)
        cds_b = CDS(100, 200)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 0

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)

    def test_len_no_overlap_right(self):
        """Tests the len_overlap function for two non-overlapping regions

        a > b"""
        cds_a = CDS(200, 300)
        cds_b = CDS(100, 200)

        hsp_a = HSP(cds_a, "", 0)
        hsp_b = HSP(cds_b, "", 0)

        expected_result = 0

        actual_result = len_overlap(hsp_a, hsp_b)

        self.assertEqual(expected_result, actual_result)
