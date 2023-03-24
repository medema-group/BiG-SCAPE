"""Contains tests for the BiG-SCAPE Run class and functions"""

# from python
from unittest import TestCase

# from dependencies
# from other modules
from src.parameters import Run, cmd_parser, Output
from src.errors import InvalidInputArgError


class TestRun(TestCase):
    """Test class for Run parameters parsing tests"""

    def setUp(self) -> None:
        self.parser = cmd_parser.cmd_parser()

    def test_parser(self):
        """Tests whether the parser object is built correctly"""

        parsed = self.parser.parse_args(
            [
                "--label",
                "test",
                "--inputdir",
                "/test/in",
                "--outputdir",
                "/test/out",
                "--pfam_dir",
                "test/pfam",
            ]
        )

        self.assertEqual(parsed.label, "test")

    def test_build_run(self):
        """Tests whether the BiG-SCAPE run object is instatiated correclty"""

        parsed = self.parser.parse_args(
            ["--inputdir", ".", "--outputdir", ".", "--pfam_dir", "."]
        )

        run = Run()
        run.parse(parsed)

        self.assertIsInstance(run, Run)

    def test_build_a_run_child(self):
        """Tests wether a run obj child is instatiated correctly"""

        parsed = self.parser.parse_args(
            ["--inputdir", ".", "--outputdir", ".", "--pfam_dir", "."]
        )
        run = Run()
        run.parse(parsed)
        output = run.output

        self.assertIsInstance(output, Output)

    def test_get_gcf_cutoffs(self):
        """Tests whether gcf cutoffs are extracted correclty"""

        parsed = self.parser.parse_args(
            [
                "--inputdir",
                ".",
                "--outputdir",
                ".",
                "--pfam_dir",
                ".",
                "--gcf_cutoffs",
                "0.1,0.2",
            ]
        )
        run = Run()
        run.parse(parsed)
        expected_cutoffs = [0.1, 0.2]

        self.assertEqual(expected_cutoffs, run.networking.gcf_cutoffs)

    def test_wrong_type_gcf_cutoffs(self):
        """Tests whether an error is raised with wrong type gcf cutoffs"""

        parsed = self.parser.parse_args(
            [
                "--inputdir",
                ".",
                "--outputdir",
                ".",
                "--pfam_dir",
                ".",
                "--gcf_cutoffs",
                "a,0.2",
            ]
        )
        run = Run()

        self.assertRaises(InvalidInputArgError, run.parse, parsed)

    def test_include_gbk(self):
        """Tests whether include gbk strings are parsed correclty"""

        parsed = self.parser.parse_args(
            [
                "--inputdir",
                ".",
                "--outputdir",
                ".",
                "--pfam_dir",
                ".",
                "--include_gbk",
                "1",
                "b",
            ]
        )
        run = Run()
        run.parse(parsed)
        expected_includes = ["1", "b"]

        self.assertEqual(expected_includes, run.input.include_gbk)
