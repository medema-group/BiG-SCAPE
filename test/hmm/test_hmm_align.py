"""Contains tests for the alignment step using pyhmmer"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.data import DB
from big_scape.genbank import GBK, CDS
from big_scape.hmm import HMMer, HSP, HSPAlignment

import big_scape.enums as bs_enums


class TestHMMAlign(TestCase):
    """Contains test methods of hmm alignment"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

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
        cds.orf_num = 1

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        hsp = HSP(cds, "PF00457.19", 249.32315063476562, 0, len(aa_seq))

        HMMer.init(hmm_path, False)

        expected_hsp = hsp
        expected_alignment = (
            "-DGMYYSFWTDGGGSVSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPL"
            "VEYYIVDNWGSYRP--TGTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGTITTGNHFDA"
            "WARAGMNMGQFRYMIMATEGYQSSGSSNIT"
        )

        expected_result = HSPAlignment(expected_hsp, expected_alignment)

        HMMer.align_simple([hsp])

        actual_result = hsp.alignment

        self.assertEqual(expected_result, actual_result)

    def test_save_hmmalignment(self):
        """Tests whether an HSP Alignment can be correctly saved to the database"""
        DB.create_in_mem()

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        gbk = GBK("", bs_enums.SOURCE_TYPE.QUERY)
        gbk.metadata = {"organism": "test", "taxonomy": "test", "description": "test"}
        gbk.save()

        cds = CDS(0, len(aa_seq) * 3)
        cds.strand = 1
        cds.orf_num = 1
        cds.aa_seq = aa_seq
        cds.parent_gbk = gbk
        cds.save()

        hsp = HSP(cds, "PF00457.19", 249.32315063476562, 0, 0)
        hsp.save()

        hsp_alignment = HSPAlignment(hsp, "")
        hsp_alignment.save()

        expected_row_count = 1

        actual_row_count = 0

        cursor_result = DB.execute_raw_query("SELECT * FROM hsp_alignment;")
        actual_row_count += len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)
