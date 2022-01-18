from unittest import TestCase

from src import big_scape

class MockOptions:
    """Test options class to replace the generated options class from argument parsing
    """
    label: str

class TestRunBase(TestCase):
    """Class containing tests for the run object, which tracks run details and keeps timing
    """
    test_run: big_scape.Run

    def setUp(self):
        self.test_run = big_scape.Run()
        self.test_run.run_mode = "Test"

        # test options object
        self.test_run.options = MockOptions()
        self.test_run.options.label = "Test"

    def test_start(self):
        """Tests whether after starting, a starting time exists in the run object
        """
        self.test_run.start(True)

        self.assertIsNotNone(self.test_run.start_time)

    def test_end(self):
        """Tests whether after starting and ending, an end time exists in run_data
        """
        self.test_run.start(True)
        self.test_run.end()

        self.assertIsNotNone(self.test_run.run_data["end_time"])
