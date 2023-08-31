"""Contains tests that check if the partial data state functions correctly determine and
return the state of partial analyses so that they can be continued
"""

# from python
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
from src.hmm import HSP, HSPAlignment


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

        actual_min_task = find_minimum_task()

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_scan(self):
        self.skipTest("Not implemented")
        expected_min_task = TASK.HMM_ALIGN

        actual_min_task = find_minimum_task()

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_align(self):
        self.skipTest("Not implemented")
        expected_min_task = TASK.HMM_ALIGN

        actual_min_task = find_minimum_task()

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_comparison(self):
        self.skipTest("Not implemented")
        expected_min_task = TASK.COMPARISON

        actual_min_task = find_minimum_task()

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
        expected_state = INPUT_TASK.NO_DATA

        actual_state = get_input_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_partial_data(self):
        self.skipTest("Not implemented")
        gbks = []
        gbk = GBK("", "test")
        gbk.add_cds_overlap_filter(CDS(0, 100))
        gbks.append(gbk)

        expected_state = INPUT_TASK.PARTIAL_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_data(self):
        self.skipTest("Not implemented")
        expected_state = INPUT_TASK.NEW_DATA

        actual_state = get_input_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_no_new_data(self):
        self.skipTest("Not implemented")
        expected_state = INPUT_TASK.SAME_DATA

        actual_state = get_input_data_state()

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
        expected_state = HMM_TASK.NO_DATA

        actual_state = get_hmm_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_new_scans_to_do(self):
        self.skipTest("Not implemented")
        expected_state = HMM_TASK.NEED_SCAN

        actual_state = get_hmm_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_no_new_scans(self):
        self.skipTest("Not implemented")
        expected_state = HMM_TASK.ALL_SCANNED

        actual_state = get_hmm_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_new_align_to_do(self):
        self.skipTest("Not implemented")
        expected_state = HMM_TASK.NEED_ALIGN

        actual_state = get_hmm_data_state()

        self.assertEqual(expected_state, actual_state)

    def test_no_new_align(self):
        self.skipTest("Not implemented")
        expected_state = HMM_TASK.ALL_ALIGNED

        actual_state = get_hmm_data_state()

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
        expected_state = COMPARISON_TASK.NO_DATA

        actual_state = get_comparison_data_state()

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
