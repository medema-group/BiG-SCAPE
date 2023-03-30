"""Contains tests for the alignment step using pyhmmer"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from src.genbank import CDS
from src.hmm import HMMer, HSP, HSPAlignment


class TestHMMAlign(TestCase):
    """Contains test methods of hmm alignment"""

    def test_align(self):
        """Tests whether a single alignment is correctly performed on a hsp"""
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

        hsp = HSP(cds, "PF00457.19", 249.32315063476562)

        HMMer.init(hmm_path, False)

        expected_hsp = hsp
        expected_alignment = (
            "-DGMYYSFWTDGGGSVSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPL"
            "VEYYIVDNWGSYRP--TGTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGTITTGNHFDA"
            "WARAGMNMGQFRYMIMATEGYQSSGSSNIT"
        )

        expected_result = HSPAlignment(expected_hsp, expected_alignment)

        actual_result = next(HMMer.align_simple([hsp]))[0]

        self.assertEqual(expected_result, actual_result)
