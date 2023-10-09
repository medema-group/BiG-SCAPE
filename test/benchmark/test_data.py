"""Tests for benchmark data loading and storage"""

# from python
from pathlib import Path
from unittest import TestCase

# from other modules
from benchmark.data import BenchmarkData


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

    def test_computed_gcf_loading(self):
        """Test loading of computed GCF assignments form database"""
        db_path = Path("test/test_data/database/valid_family.db")

        expected_data = {
            0.3: {
                "Chr1.cluster024": 1,
                "CM000574.1.cluster001": 2,
                "CM000574.1.cluster030": 3,
                "CM000578.1.cluster010": 11,
                "CM000602.2.cluster010": 11,
                "CM003198.1.cluster026": 11,
                "HF679024.1.cluster002": 11,
                "HG323944.1.cluster001": 2,
                "HG323944.1.cluster030": 11,
                "JH717896.1.cluster018": 11,
            }
        }

        dataloader = BenchmarkData(None, db_path)
        dataloader.load_computed_labels()

        self.assertEqual(dataloader.computed_labels, expected_data)
