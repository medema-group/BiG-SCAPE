"""module containing integration tests for loading files and records"""

# from python
import os
from pathlib import Path
from unittest import TestCase

# from other modules
import big_scape.enums as bs_enums
import big_scape.data as bs_data
import big_scape.hmm as bs_hmm
from big_scape.genbank import GBK, CDS

# from this module
from big_scape.file_input.load_files import (
    load_gbks,
    load_dataset_folder,
    remove_duplicate_gbk,
    bgc_length_contraint,
)


class TestLoadRecords(TestCase):
    """Test class for loading files and records in files"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_load_files_empty_db(self):
        """Tests whether the load_dataset_folder function correctly loads files into an empty database"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "query_bgc_path": None,
            "mibig_version": None,
            "reference_dir": None,
            "db_path": Path("test/test_data/database/not_exists.db"),
            "include_gbk": ["*"],
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        if not run["db_path"].exists():
            bs_data.DB.create_in_mem()

        bigscape_dir = Path(os.path.dirname(os.path.abspath(__file__)))

        input_gbks = load_gbks(run, bigscape_dir)

        self.assertIsNotNone(input_gbks)

    def test_load_gbks_load_dataset_simple(self):
        """Tests whether the simplest workflow of loading gbks from a folder and performing filter/checks works"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "query_bgc_path": None,
            "mibig_version": None,
            "reference_dir": None,
            "db_path": Path("test/test_data/database/not_exists.db"),
            "include_gbk": ["gbk"],
            "exclude_gbk": ["final"],
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
        }

        input_gbks = []

        gbks = load_dataset_folder(run["input_dir"], run, bs_enums.SOURCE_TYPE.QUERY)
        input_gbks.extend(gbks)

        self.assertIsNotNone(input_gbks)

        length = len(input_gbks)

        # we have 4 gbks in the test folder and one of them is named "final", thus excluded
        self.assertEqual(length, 3)

        input_gbks = remove_duplicate_gbk(input_gbks)

        length = len(input_gbks)

        # we have 3 gbks in input but they are all the same
        self.assertEqual(length, 1)

        input_gbks = bgc_length_contraint(input_gbks)

        length = len(input_gbks)

        # that one gbk passes the length constraint
        self.assertEqual(length, 1)

    def test_load_gbks_load_dataset_reference(self):
        """Tests whether the simplest workflow of loading gbks works"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "query_bgc_path": None,
            "mibig_version": None,
            "reference_dir": Path("test/test_data/valid_gbk_folder/"),
            "db_path": Path("test/test_data/database/not_exists.db"),
            "include_gbk": ["gbk"],
            "exclude_gbk": ["final"],
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
            "force_gbk": False,
        }

        input_gbks = []
        gbks_query = load_dataset_folder(
            run["input_dir"], run, bs_enums.SOURCE_TYPE.QUERY
        )
        input_gbks.extend(gbks_query)

        gbks_ref = load_dataset_folder(
            run["reference_dir"], run, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(gbks_ref)

        self.assertIsNotNone(input_gbks)

        length = len(input_gbks)

        # we have 4 gbks in the input folder but one is named "final", thus excluded
        # we have 4 gbks in the reference folder
        self.assertEqual(length, 7)

        input_gbks = remove_duplicate_gbk(input_gbks)
        input_gbks = bgc_length_contraint(input_gbks)

        length = len(input_gbks)

        # in reality there are 4 unique gbks, and they all pass the length constraint
        self.assertEqual(length, 4)

        has_refs = False
        has_queries = False

        for gbk in input_gbks:
            if gbk.source_type == bs_enums.SOURCE_TYPE.REFERENCE:
                has_refs = True
            if gbk.source_type == bs_enums.SOURCE_TYPE.QUERY:
                has_queries = True

        self.assertTrue(has_refs)
        self.assertTrue(has_queries)

    def test_load_gbks_find_minimum_task(self):
        """Tests whether find_minimum_task correctly identifies the task to do next"""

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "query_bgc_path": None,
            "mibig_version": None,
            "reference_dir": None,
            "db_path": Path("test/test_data/database/not_exists.db"),
            "include_gbk": ["gbk"],
            "exclude_gbk": ["final"],
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
            "force_gbk": False,
        }

        if not run["db_path"].exists():
            bs_data.DB.create_in_mem()

        input_gbks = []
        gbks_query = load_dataset_folder(
            run["input_dir"], run, bs_enums.SOURCE_TYPE.QUERY
        )
        input_gbks.extend(gbks_query)

        input_gbks = remove_duplicate_gbk(input_gbks)

        task_state = bs_data.find_minimum_task(input_gbks)

        self.assertEqual(task_state, bs_enums.TASK.SAVE_GBKS)

        missing_gbks = bs_data.get_missing_gbks(input_gbks)

        self.assertEqual(len(missing_gbks), 4)

        gbk1 = missing_gbks[0]
        gbk1.save_all()

        task_state = bs_data.find_minimum_task(input_gbks)

        self.assertEqual(task_state, bs_enums.TASK.SAVE_GBKS)

        gbk2 = missing_gbks[1]
        gbk2.save_all()
        gbk3 = missing_gbks[2]
        gbk3.save_all()
        gbk4 = missing_gbks[3]
        gbk4.save_all()

        task_state = bs_data.find_minimum_task(input_gbks)

        self.assertNotEqual(task_state, bs_enums.TASK.SAVE_GBKS)
        self.assertEqual(task_state, bs_enums.TASK.HMM_SCAN)

    def test_load_gbks_save_to_and_load_from_db(self):
        """Tests whether gbks are correclty saved to db, and then loaded from db"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "query_bgc_path": None,
            "mibig_version": None,
            "reference_dir": Path("test/test_data/valid_gbk_folder/"),
            "db_path": Path("test/test_data/database/not_exists.db"),
            "include_gbk": ["gbk"],
            "exclude_gbk": ["final"],
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "legacy_classify": False,
            "force_gbk": False,
        }

        if not run["db_path"].exists():
            bs_data.DB.create_in_mem()

        input_gbks = []
        gbks_query = load_dataset_folder(
            run["input_dir"], run, bs_enums.SOURCE_TYPE.QUERY
        )
        input_gbks.extend(gbks_query)

        gbks_ref = load_dataset_folder(
            run["reference_dir"], run, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(gbks_ref)

        input_gbks = remove_duplicate_gbk(input_gbks)

        missing_gbks = bs_data.get_missing_gbks(input_gbks)

        for gbk in missing_gbks:
            gbk.save_all()

        task_state = bs_data.find_minimum_task(input_gbks)
        self.assertNotEqual(task_state, bs_enums.TASK.SAVE_GBKS)

        source_dict = {gbk.hash: gbk.source_type for gbk in input_gbks}

        input_gbks_from_db = GBK.load_many(input_gbks)
        for gbk in input_gbks_from_db:
            gbk.source_type = source_dict[gbk.hash]
            bs_hmm.HSP.load_all(gbk.genes)

        self.assertEqual(len(input_gbks), len(input_gbks_from_db))

        gbk1 = input_gbks_from_db[0]

        self.assertIsInstance(gbk1, GBK)
        self.assertIsInstance(gbk1.genes[0], CDS)

        for gbk in input_gbks:
            self.assertIn(gbk, input_gbks_from_db)

            idx = input_gbks_from_db.index(gbk)
            gbk_db = input_gbks_from_db[idx]
            self.assertEqual(gbk_db.source_type, gbk.source_type)
