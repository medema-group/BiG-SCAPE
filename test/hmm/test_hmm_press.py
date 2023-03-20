"""Contains tests for pressing hmm files into other formats used in scanning and aligning"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies

# from other modules
from src.hmm import hmm_press


class TestHMMPress(TestCase):
    """Contains tests to make sure hmm files can be correctly pressed into optimized formats"""

    def test_press_input(self):
        """Tests whether the pressing process completes successfully"""
        input_hmm_path = Path("test/test_data/hmm/valid.hmm")

        hmm_press(input_hmm_path)

        expected_file_paths = [
            Path("test/test_data/hmm/valid.hmm.h3f"),
            Path("test/test_data/hmm/valid.hmm.h3i"),
            Path("test/test_data/hmm/valid.hmm.h3m"),
            Path("test/test_data/hmm/valid.hmm.h3p"),
        ]

        existing_paths = [path.exists for path in expected_file_paths]

        self.assertTrue(all(existing_paths))

    def test_press_input_invalid(self):
        """Tests whether press_hmm raises a ValueError when an invalid hmm file is given"""
        input_hmm_path = Path("test/test_data/hmm/invalid.hmm")

        self.assertRaises(ValueError, hmm_press, input_hmm_path)
