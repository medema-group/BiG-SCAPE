"""Contains tests for setting up and breaking down the HMMer class variables"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from pyhmmer.plan7 import HMM, OptimizedProfile

# from other modules
from src.hmm import HMMer


class TestHMMSetup(TestCase):
    """Contains tests for setting up and breaking down the HMMer class variables"""

    def class_cleanup(sefl):
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

        profile = HMMer.get_profile("PF00457.19")

        self.assertIsInstance(profile, OptimizedProfile)
