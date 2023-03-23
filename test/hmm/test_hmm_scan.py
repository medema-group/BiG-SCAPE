"""Contains tests to test the hmm scanning functionality using pyhmmer"""

# from python
from pathlib import Path
from unittest import TestCase

# from other modules
from src.hmm import HMMer
from src.genbank import CDS
from src.hmm import HSP


class TestHMMScan(TestCase):
    """Contains tests to check the hmmscan functionality"""

    def test_scan_worker_method(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""
        pass
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

        actual_result = next(HMMer.scan([cds]))

        self.assertEqual(expected_result, actual_result)
