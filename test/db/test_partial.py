"""Contains tests that check if the partial data state functions correctly determine and
return the state of partial analyses so that they can be continued
"""

# from python
from pathlib import Path
import random
from unittest import TestCase

# from other modules
from src.data import (
    DB,
    find_minimum_task,
    get_input_data_state,
    get_hmm_data_state,
    get_comparison_data_state,
    TASK,
    INPUT_TASK,
    HMM_TASK,
    COMPARISON_TASK,
)
from src.genbank import GBK, CDS
from src.hmm import HSP, HSPAlignment, HMMer


def create_mock_gbk() -> GBK:
    gbk = GBK(Path(f"test_path_{random.randint(0, 1000)}.gbk"), "test")
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    return gbk


def add_mock_hsp_cds(cds: CDS) -> None:
    hsp = HSP(cds, "PF01234.12", 1.0, 0, 100)
    cds.hsps.append(hsp)


def add_mock_hsp_alignment_hsp(hsp: HSP) -> None:
    hsp_alignment = HSPAlignment(hsp, "AAAAAAAAAA")
    hsp.alignment = hsp_alignment


class TestPartial(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_min_task_data(self):
        DB.create_in_mem()
        expected_min_task = TASK.LOAD_GBKS

        gbks = [create_mock_gbk()]

        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_scan(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]
        # add gbk to db
        gbks[0].save_all()

        expected_min_task = TASK.HMM_SCAN
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_align(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]
        # add hsp to gbk
        add_mock_hsp_cds(gbks[0].genes[0])

        # add gbk to db
        gbks[0].save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])

        expected_min_task = TASK.HMM_ALIGN
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_comparison(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]
        # add hsp to gbk
        add_mock_hsp_cds(gbks[0].genes[0])
        # add hspalignment to hsp
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])

        # add gbk to db
        gbks[0].save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])
        # sigh
        gbks[0].genes[0].hsps[0].alignment.save()

        expected_min_task = TASK.COMPARISON
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)


class TestPartialInputs(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_data(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]

        expected_state = INPUT_TASK.NO_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_partial_data(self):
        DB.create_in_mem()

        # gbk 1
        gbk_1_in_db = create_mock_gbk()
        gbk_1_in_db.save_all()
        # gbk 2
        gbk_2_in_db = create_mock_gbk()
        gbk_2_in_db.save_all()

        gbks = [gbk_1_in_db]

        expected_state = INPUT_TASK.PARTIAL_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_data(self):
        DB.create_in_mem()

        gbk_in_db = create_mock_gbk()
        gbk_in_db.save_all()

        gbk_not_in_db = create_mock_gbk()
        gbks = [gbk_in_db, gbk_not_in_db]

        expected_state = INPUT_TASK.NEW_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_no_new_data(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(),
            create_mock_gbk(),
            create_mock_gbk(),
        ]
        for gbk in gbks:
            gbk.save_all()

        expected_state = INPUT_TASK.SAME_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)


class TestPartialHMM(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_scans_done(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]

        expected_state = HMM_TASK.NO_DATA
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_scans_to_do(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(),
            create_mock_gbk(),
        ]
        # add hsp to one gbk
        add_mock_hsp_cds(gbks[0].genes[0])

        # add gbks to db
        for gbk in gbks:
            gbk.save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])

        expected_state = HMM_TASK.NEED_SCAN
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_align_to_do(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(),
            create_mock_gbk(),
        ]
        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        # add hspalignment to only one hsp
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])

        # add gbks and hsps to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
        # add hspalignment to db
        gbks[0].genes[0].hsps[0].alignment.save()

        expected_state = HMM_TASK.NEED_ALIGN
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_no_new_align(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(),
            create_mock_gbk(),
        ]
        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        # add hspalignment to both hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        expected_state = HMM_TASK.ALL_ALIGNED
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)


class TestPartialComparison(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_comparisons_done(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk()]

        expected_state = COMPARISON_TASK.NO_DATA
        actual_state = get_comparison_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_comparisons_to_do(self):
        self.skipTest("Not implemented")
        expected_state = COMPARISON_TASK.NEW_DATA

        actual_state = get_comparison_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_no_new_comparisons(self):
        self.skipTest("Not implemented")
        expected_state = COMPARISON_TASK.ALL_DONE

        actual_state = get_comparison_data_state()

        self.assertEqual(expected_state, actual_state)
