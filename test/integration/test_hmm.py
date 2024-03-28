"""module containing integration tests for the hmm search/scan and align workflow"""

# from python
from pathlib import Path
from unittest import TestCase
import platform
import unittest

# from other modules
import big_scape.data as bs_data
import big_scape.hmm as bs_hmm
from big_scape.genbank import GBK, CDS, Region
import big_scape.enums as bs_enums


def create_mock_gbk(i) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), bs_enums.SOURCE_TYPE.QUERY)
    aa_seq = (
        "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
        "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
        "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
        "IMATEGYQSSGSSNITVSG"
    )
    cds = CDS(0, len(aa_seq) * 3)
    cds.aa_seq = aa_seq
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, 1, 0, 100, False, "test")
    gbk.region._db_id = i
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def add_mock_hsp_cds(cds: CDS) -> None:
    hsp = bs_hmm.HSP(cds, "PF01234", 1.0, 0, 100)
    cds.hsps.append(hsp)


def add_mock_hsp_alignment_hsp(hsp: bs_hmm.HSP) -> None:
    hsp_alignment = bs_hmm.HSPAlignment(hsp, "AAAAAAAAAA")
    hsp.alignment = hsp_alignment


class TestHMMScan(TestCase):
    """Test class for hmm search/scan workflow"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def class_cleanup(self):
        bs_hmm.profiles = []
        bs_hmm.profile_index = {}
        bs_hmm.pipeline = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)
        self.addCleanup(self.class_cleanup)

    def test_hmmsearch_simple_workflow(self):
        """Tests the entire hmmsearch worfkflow"""

        bs_data.DB.create_in_mem()

        gbk = create_mock_gbk(1)
        gbk.save_all()

        cds_to_scan = bs_data.get_cds_to_scan([gbk])
        expected_cds_to_scan = [gbk.genes[0]]

        self.assertEqual(expected_cds_to_scan, cds_to_scan)

        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")
        bs_hmm.HMMer.init(hmm_path)

        bs_hmm.HMMer.hmmsearch_simple(cds_to_scan, cores=1)

        expected_result = bs_hmm.HSP(
            cds_to_scan[0], "PF00457", 249.32315063476562, 0, 0
        )
        actual_result = cds_to_scan[0].hsps[0]

        self.assertEqual(expected_result, actual_result)

        for cds in cds_to_scan:
            for hsp in cds.hsps:
                hsp.save(False)

        cds_to_scan = bs_data.get_cds_to_scan([gbk])
        expected_cds_to_scan = []

        self.assertEqual(expected_cds_to_scan, cds_to_scan)

        run_state = bs_data.find_minimum_task([gbk])
        self.assertEqual(run_state, bs_enums.TASK.HMM_ALIGN)

    @unittest.skipIf(platform.system() == "Darwin", "Unsupported OS")
    def test_hmmsearch_multiprocess_workflow(self):
        """Tests the entire hmmsearch worfkflow"""

        bs_data.DB.create_in_mem()

        gbk = create_mock_gbk(1)
        gbk.save_all()

        cds_to_scan = bs_data.get_cds_to_scan([gbk])
        expected_cds_to_scan = [gbk.genes[0]]

        self.assertEqual(expected_cds_to_scan, cds_to_scan)

        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")
        bs_hmm.HMMer.init(hmm_path)

        bs_hmm.HMMer.hmmsearch_multiprocess(cds_to_scan, batch_size=1, cores=1)

        expected_result = bs_hmm.HSP(
            cds_to_scan[0], "PF00457", 249.32315063476562, 0, 0
        )
        actual_result = cds_to_scan[0].hsps[0]

        self.assertEqual(expected_result, actual_result)

        for cds in cds_to_scan:
            for hsp in cds.hsps:
                hsp.save(False)

        cds_to_scan = bs_data.get_cds_to_scan([gbk])
        expected_cds_to_scan = []

        self.assertEqual(expected_cds_to_scan, cds_to_scan)

        run_state = bs_data.find_minimum_task([gbk])
        self.assertEqual(run_state, bs_enums.TASK.HMM_ALIGN)

    def test_hmmscan_workflow_file(self):
        """Tests the entire hmm scan workflow with a file as input"""

        bs_data.DB.create_in_mem()

        gbk = create_mock_gbk(1)
        gbk.save_all()

        cds_to_scan = bs_data.get_cds_to_scan([gbk])

        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")
        bs_hmm.HMMer.init(hmm_path, False)

        bs_hmm.HMMer.hmmsearch_simple(cds_to_scan, cores=1)

        for cds in cds_to_scan:
            for hsp in cds.hsps:
                hsp.save(False)

        run_state = bs_data.find_minimum_task([gbk])
        self.assertEqual(run_state, bs_enums.TASK.HMM_ALIGN)

        hsps_to_align = bs_data.get_hsp_to_align([gbk])
        expected_hsps_to_align = [cds_to_scan[0].hsps[0]]

        self.assertEqual(expected_hsps_to_align, hsps_to_align)

        bs_hmm.HMMer.align_simple(hsps_to_align)

        expected_alignment = (
            "-DGMYYSFWTDGGGSVSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPL"
            "VEYYIVDNWGSYRP--TGTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGTITTGNHFDA"
            "WARAGMNMGQFRYMIMATEGYQSSGSSNIT"
        )

        expected_result = bs_hmm.HSPAlignment(hsps_to_align[0], expected_alignment)

        hsp = hsps_to_align[0]
        actual_result = hsp.alignment

        run_state = bs_data.find_minimum_task([gbk])

        self.assertEqual(expected_result, actual_result)

        for hsp in hsps_to_align:
            if hsp.alignment is None:
                continue
            hsp.alignment.save(False)

        cds_to_scan = bs_data.get_cds_to_scan([gbk])
        expected_cds_to_scan = []
        self.assertEqual(expected_cds_to_scan, cds_to_scan)

        run_state = bs_data.find_minimum_task([gbk])
        self.assertEqual(run_state, bs_enums.TASK.COMPARISON)
