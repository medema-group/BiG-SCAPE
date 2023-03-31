"""Contains tests for pressing hmm files into other formats used in scanning and aligning"""

# from python
from pathlib import Path
from unittest import TestCase


# from other modules
from src.hmm import HMMer


class TestHMMPress(TestCase):
    """Contains tests to make sure hmm files can be correctly pressed into optimized formats"""

    def remove_pressed_files(self):
        """Removes any pressed files from these tests in the test_data directory"""
        for created_file in self.created_files:
            created_file.unlink(True)

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.created_files: list[Path] = []

        self.addCleanup(self.remove_pressed_files)

    def test_press_input(self):
        """Tests whether the pressing process correctly generates pressed hmm files"""
        input_hmm_path = Path("test/test_data/hmm/valid.hmm")

        HMMer.press(input_hmm_path)

        expected_file_paths = [
            Path("test/test_data/hmm/valid.hmm.h3f"),
            Path("test/test_data/hmm/valid.hmm.h3i"),
            Path("test/test_data/hmm/valid.hmm.h3m"),
            Path("test/test_data/hmm/valid.hmm.h3p"),
        ]

        self.created_files.extend(expected_file_paths)

        existing_paths = [path.exists for path in expected_file_paths]

        self.assertTrue(all(existing_paths))

    def test_press_input_invalid(self):
        """Tests whether press_hmm raises a ValueError when an invalid hmm file is given"""
        input_hmm_path = Path("test/test_data/hmm/invalid.hmm")

        self.assertRaises(ValueError, HMMer.press, input_hmm_path)
