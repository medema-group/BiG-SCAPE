"""Contains tests to test the hmm scanning functionality using pyhmmer"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from pyhmmer.easel import TextSequence

# from other modules
from src.hmm.hmmer import HMMer, cds_to_input_task
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

        expected_result = HSP(cds, "PF00457.19", 249.32315063476562, 0, 0)

        HMMer.init(hmm_path)

        actual_result = next(HMMer.hmmsearch_simple([cds]))

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

        expected_result = HSP(cds, "PF00457.19", 249.32315063476562, 0, 0)

        actual_result = list(HMMer.hmmsearch_multiprocess([cds]))[0]

        self.assertEqual(expected_result, actual_result)
