"""Contains tests for sourmash utilities."""

# from python
from unittest import TestCase
from pathlib import Path
import tempfile
import os

# from dependencies

# from other modules
import big_scape.enums as bs_enums
from big_scape.dereplicating.gbk_components.gbk import GBK
from big_scape.dereplicating.gbk_components.cds import CDS
from big_scape.dereplicating.networking import Edge
from big_scape.dereplicating.sourmash_utilities import (
    get_fasta_names_paths,
    write_concat_cds_fasta,
    make_sourmash_input,
    parse_sourmash_results,
)


class TestSourmashUtilities(TestCase):
    """Test class for sourmash utilities."""

    def test_make_sourmash_input(self):
        """Tests whether the make_sourmash_input function correctly generates sourmash input files."""
        with tempfile.TemporaryDirectory() as tmpdir:

            run = {
                "input_dir": Path("input_dir/"),
                "output_dir": Path(tmpdir),
                "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
                "include_gbk": ["cluster", "region"],
                "exclude_gbk": ["final"],
                "mode": bs_enums.RUN_MODE.DEREPLICATE,
                "cores": 1
            }

            gbk_1 = GBK(Path("input_dir/gbk_1.gbk"), "hash_1", 10, "1", bs_enums.SOURCE_TYPE.QUERY)
            cds_1 = CDS(0, 1, 1, None, "A"*2)
            cds_2 = CDS(0, 1, 1, None, "B"*2)
            gbk_1.components[bs_enums.COMPONENTS.CDS] = [cds_1, cds_2]

            gbk_list = [gbk_1]

            sourmash_dir, CDS_fasta_dir, manysketch_csv_path = make_sourmash_input(gbk_list, run)

            # Check if the directories were created
            self.assertTrue(os.path.exists(sourmash_dir))
            self.assertTrue(os.path.exists(CDS_fasta_dir))

            expected_sourmash_dir = Path(os.path.join(run["output_dir"], "sourmash"))
            expected_CDS_fasta_dir = Path(os.path.join(expected_sourmash_dir, "CDS_fastas"))

            self.assertEqual(sourmash_dir, expected_sourmash_dir)
            self.assertEqual(CDS_fasta_dir, expected_CDS_fasta_dir)

            # Check if the fasta files were created
            fasta_files = list(Path(CDS_fasta_dir).glob("*.fasta"))
            self.assertEqual(len(fasta_files), 1)
            self.assertTrue(os.path.exists(fasta_files[0]))
            expected_fasta_path = Path(os.path.join(CDS_fasta_dir, "gbk_1.gbk.hash_1.fasta"))
            self.assertEqual(fasta_files[0], expected_fasta_path)

            # Check if the manysketch cvs file was created
            expected_manysketch_csv_path = Path(os.path.join(sourmash_dir, "manysketch.csv"))
            self.assertEqual(manysketch_csv_path, expected_manysketch_csv_path)
            self.assertTrue(os.path.exists(manysketch_csv_path))

            # Check if the manysketch csv file contains the correct data
            with open(manysketch_csv_path, "r") as manysketch_csv_file:
                header = manysketch_csv_file.readline().strip()
                expected_header = "name,genome_filename,protein_filename"
                self.assertEqual(header, expected_header)

                line_1 = manysketch_csv_file.readline().strip()
                expected_line_1 = f"gbk_1.gbk.hash_1,,{expected_fasta_path}"
                self.assertEqual(line_1, expected_line_1)

    def test_making_fasta_names(self):
        """Test the making of fasta names."""

        run = {
            "input_dir": Path("data_input_dir/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["cluster", "region"],
            "exclude_gbk": ["final"],
            "mode": bs_enums.RUN_MODE.DEREPLICATE,
            "cores": 1
        }

        gbk_1 = GBK(Path("data_input_dir/gbk_folder/gbk_1.gbk"), "gbk_hash", 10, "as_version", bs_enums.SOURCE_TYPE.QUERY)

        sourmash_dir = Path("soumarsh_folder/sourmash_input/")

        gbk_name, fasta_path = get_fasta_names_paths(gbk_1, run, sourmash_dir)

        expected_fasta_name = Path("soumarsh_folder/sourmash_input/gbk_folder.gbk_1.gbk.gbk_hash.fasta")
        expected_gbk_name = "gbk_folder.gbk_1.gbk"

        self.assertEqual(fasta_path, expected_fasta_name)
        self.assertEqual(gbk_name, expected_gbk_name)

    def test_write_concat_fasta(self):
        """Test the writing of concatenated fasta files."""

        with tempfile.TemporaryDirectory() as tmpdir:
            sourmash_dir = tmpdir

            run = {
                "input_dir": Path("input_dir/"),
                "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
                "include_gbk": ["cluster", "region"],
                "exclude_gbk": ["final"],
                "mode": bs_enums.RUN_MODE.DEREPLICATE,
                "cores": 1
            }

            cds_1 = CDS(0, 1, 1, None, "A"*2)
            cds_2 = CDS(0, 1, 1, None, "B"*2)

            gbk = GBK(Path("input_dir/gbk_1.gbk"), "hash", 10, "1", bs_enums.SOURCE_TYPE.QUERY)
            gbk.components[bs_enums.COMPONENTS.CDS] = [cds_1, cds_2]

            CDS.concatenate_cds(gbk)

            gbk_name, fasta_path = get_fasta_names_paths(gbk, run, sourmash_dir)

            write_concat_cds_fasta(gbk, gbk_name, fasta_path)

            # Check if the fasta file was created
            self.assertTrue(os.path.exists(fasta_path))

            # Check if the contents of the fasta file are correct
            with open(fasta_path, "r") as fasta_file:
                header = fasta_file.readline().strip()
                seq = fasta_file.readline().strip()
                expected_header = f">{gbk_name}.{gbk.hash}"
                expected_seq = "A"*2 + "B"*2
                self.assertEqual(header, expected_header)
                self.assertEqual(seq, expected_seq)

    def test_parse_sourmash_results(self):
        """Test the parsing of sourmash results."""

        sourmash_pairwise_csv = Path("test/test_data/sourmash_output/sourmash_branchwater_pairwise_output.csv")

        edges = parse_sourmash_results(sourmash_pairwise_csv, 0.3)

        self.assertEqual(len(edges), 3)

        expected_edges = [
            Edge("gbk_1", "gbk_2", 0.31),
            Edge("gbk_3", "gbk_2", 0.35),
            Edge("gbk_1", "gbk_3", 0.70)
        ]

        self.assertIn(expected_edges[0], edges)
        self.assertIn(expected_edges[1], edges)
        self.assertIn(expected_edges[2], edges)

