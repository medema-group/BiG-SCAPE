"""module containing integration tests for the family calling workflows"""

# from python
from unittest import TestCase

# from dependencies

# from other modules
import big_scape.data as bs_data

# from this module


class TestComparison(TestCase):
    """Test class for hmm search/scan workflow"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_family_workflow_mix_one_cutoff(self):
        self.skipTest("Not implemented yet")

    def test_family_workflow_mix_two_cutoffs(self):
        self.skipTest("Not implemented yet")
