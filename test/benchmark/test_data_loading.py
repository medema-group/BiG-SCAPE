"""Tests for benchmark data loading and storage"""

# from python
import logging
from pathlib import Path
from unittest import TestCase

# from other modules
from big_scape.benchmarking import BenchmarkData


class TestBenchmarkData(TestCase):
    """Test class for data loading tasks"""

    def test_curated_gcf_loading(self):
        """Test the loading of curated GCF assignments"""
        expected_data = {
            "AmC": "II",
            "Chr1.cluster049": "Fusarium_117",
            "CM000578.1.cluster041": "Fusarium_60",
            "NC_012490.1.cluster047": "nrps_5",
            "AC-40.region18": "JK1_GCF_25",
            "BGC0000315": "JK1_GCF_26",
            "BGC0002412": "Pyrrothine antibiotics",
            "BGC0002415": "NRP siderophores",
        }
        curated_path = Path("test/test_data/curated_gcfs/valid_curated_gcfs.tsv")
        testloader = BenchmarkData(curated_path, None)
        testloader.load_curated_labels()

        self.assertEqual(testloader.curated_labels, expected_data)

    def test_empty_assignments_warning(self):
        """Test warning upon loading empty GCF assignment file"""
        empty_file = Path("test/test_data/curated_gcfs/empty_curated_gcfs.tsv")
        dataloader = BenchmarkData(None, None)
        with self.assertLogs(level=logging.INFO) as cm:
            logging.info("nonsense")
            dataloader.read_gcf_tsv(empty_file)

        # cm.output a list of strings of all the logs
        str = "GCF assignment file is empty"
        warning = any(str in log for log in cm.output)

        self.assertEqual(warning, True)

    def test_bs1_computed_gcf_loading(self):
        """Test loading of BS1 computed GCF assignments from files"""
        run_path = Path("test/test_data/bs1_output/valid_network_files")

        expected_data = {
            "0.30": {
                "BDCX01000013.1.region001": "0",
                "c00002_PLM3_2_...region001": "0",
                "BGC0000611": "2",
                "BGC0000614": "2",
            },
            "0.60": {
                "BDCX01000013.1.region001": "8",
                "BGC0000611": "8",
                "BGC0000614": "8",
                "c00002_PLM3_2_...region001": "8",
            },
        }

        dataloader = BenchmarkData(None, None)
        dataloader.load_computed_bs1_labels(run_path)
        self.assertEqual(dataloader.computed_labels, expected_data)

    def test_bs2_computed_gcf_loading(self):
        """Test loading of BS2 computed GCF assignments from database"""
        db_path = Path("test/test_data/database/valid_family.db")

        expected_data = {
            "0.3": {
                "Chr1.cluster024": "1",
                "CM000574.1.cluster001": "2",
                "CM000574.1.cluster030": "3",
                "CM000578.1.cluster010": "11",
                "CM000602.2.cluster010": "11",
                "CM003198.1.cluster001": "6",
                "CM003198.1.cluster026": "11",
                "HF679024.1.cluster002": "11",
                "HG323944.1.cluster001": "2",
                "HG323944.1.cluster030": "11",
                "JH717896.1.cluster018": "11",
            }
        }

        dataloader = BenchmarkData(None, None)
        dataloader.load_computed_bs2_labels(db_path)

        self.assertEqual(dataloader.computed_labels, expected_data)

    def test_bslice_computed_gcf_loading(self):
        """Test loading of BiG-SLICE computed GCF assignments from database"""
        db_path = Path("test/test_data/database/valid_bigslice.db")

        expected_data = {
            "0.4": {"NC_003888.3.region002": "1", "NC_003888.3.region001": "2"},
            "0.3": {"NC_003888.3.region002": "3", "NC_003888.3.region001": "4"},
            "0.6": {"NC_003888.3.region002": "5", "NC_003888.3.region001": "6"},
        }

        dataloader = BenchmarkData(None, None)
        dataloader.load_computed_bslice_labels(db_path)

        self.assertEqual(dataloader.computed_labels, expected_data)

    def test_missing_bs1_data_raises(self):
        """Test if correct errors are raised upon missing BS1 results"""
        no_runs = Path("test/test_data/bs1_output/empty_network_files")
        no_mix_data = Path("test/test_data/bs1_output/no_mix_network_files")
        dataloader = BenchmarkData(None, None)
        self.assertRaises(
            FileNotFoundError,
            dataloader.load_computed_bs1_labels,
            no_runs,
        )
        self.assertRaises(
            FileNotFoundError,
            dataloader.load_computed_bs1_labels,
            no_mix_data,
        )

    def test_missing_database_raises(self):
        """Test if correct errors are raised upon missing BS2 results"""
        no_db_path = Path("test/test_data/database")
        dataloader = BenchmarkData(None, no_db_path)
        self.assertRaises(FileNotFoundError, dataloader.load_computed_labels)