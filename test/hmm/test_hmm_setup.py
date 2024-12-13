"""Contains tests for setting up and breaking down the HMMer class variables"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from pyhmmer.plan7 import HMM, OptimizedProfile

# from other modules
from big_scape.hmm import HMMer
import big_scape.hmm as bs_hmm


class TestHMMSetup(TestCase):
    """Contains tests for setting up and breaking down the HMMer class variables"""

    def class_cleanup(self):
        HMMer.profiles = []
        HMMer.profile_index = {}
        HMMer.pipeline = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.class_cleanup)

    def test_init_optimized(self):
        """Tests whether init() correctly loads the HMMer.profiles object with optimized
        profiles
        """
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        self.assertIsInstance(HMMer.profiles[0], OptimizedProfile)

    def test_init_unoptimized(self):
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path, False)

        self.assertIsInstance(HMMer.profiles[0], HMM)

    def test_unload(self):
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        HMMer.unload()

        expected_true_checks = [
            HMMer.profiles == [],
            HMMer.profile_index == {},
            HMMer.pipeline is None,
        ]

        self.assertTrue(all(expected_true_checks))

    def test_get_profile(self):
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        profile = HMMer.get_profile("PF00457")

        self.assertIsInstance(profile, OptimizedProfile)

    def test_get_cds_batches(self):
        """tests whether cds_batch_generator correctly splits all cds in bacthes"""

        cds_list = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]
        cores = 3

        expected_batches = [["a", "b", "c"], ["d", "e", "f"], ["g", "h", "i"], ["j"]]

        actual_batches = list(bs_hmm.cds_batch_generator(cds_list, cores))

        self.assertEqual(expected_batches, actual_batches)

        cores = 1

        expected_batches = [["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]]
        actual_batches = list(bs_hmm.cds_batch_generator(cds_list, cores))

        self.assertEqual(expected_batches, actual_batches)

        cores = 5

        expected_batches = [["a", "b", "c", "d", "e", "f", "g", "h", "i", "j"]]
        actual_batches = list(
            bs_hmm.cds_batch_generator(cds_list, cores, spread_input=False)
        )

        self.assertEqual(expected_batches, actual_batches)
