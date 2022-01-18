from unittest import TestCase

from src import big_scape

class test_options:
    """Test options class to replace the generated options class from argument parsing
    """
    label: str
    inputdir: str
    outputdir: str
    pfam_dir: str
    cores: int
    include_gbk_str: str
    exclude_gbk_str: str
    verbose: bool
    include_singletons: bool
    domain_overlap_cutoff: float
    min_bgc_size: int
    mix: bool
    no_classify: bool
    banned_classes: str
    cutoffs: float
    clans: bool
    clan_cutoff: float
    hybrids: bool
    mode: str
    anchorfile: str
    force_hmmscan: bool
    skip_ma: bool
    mibig21: bool
    mibig14: bool
    mibig13: bool
    mibig_path: str
    query_bgc: str
    domain_includelist: bool


class TestRunBase(TestCase):
    """Class containing tests for the run object, which tracks run details and keeps timing
    """
    test_run: big_scape.Run

    def setUp(self):
        self.test_run = big_scape.Run()
        self.test_run.run_mode = "Test"

        # test options object
        self.test_run.options = test_options()
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
