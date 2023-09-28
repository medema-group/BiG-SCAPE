"""Contains tests for the HSP class"""


# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.genbank import GBK
from big_scape.data import DB
from big_scape.hmm import HSP


class TestHSP(TestCase):
    """Contains tests for HSP class functions"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_load_all(self):
        DB.load_from_disk(Path("test/test_data/database/valid_populated.db"))

        all_gbks = GBK.load_all()

        expected_hsp_count = 501
        actual_hsp_count = 0

        for gbk in all_gbks:
            HSP.load_all(gbk.genes)
            actual_hsp_count += sum([len(cds.hsps) for cds in gbk.genes])

        self.assertEqual(expected_hsp_count, actual_hsp_count)
