"""Tests for benchmark metrics calculation"""

# from python
from unittest import TestCase

# from other modules
from big_scape.benchmarking import BenchmarkMetrics


class TestBenchmarkMetrics(TestCase):
    """Test class for metrics calculation tasks"""

    def test_vmeasure(self):
        """Test the calculation of homogeneity, completeness, v_measure"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        hcv = calc.calculate_v_measure()
        self.assertEqual(hcv, (1, 1, 1))

    def test_purity(self):
        """Test the calculation of purity"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        p = calc.calculate_purity()

        expected = {"gcf1": 1.0, "gcf2": 1.0}
        self.assertEqual(p, expected)

    def test_entropy(self):
        """Test the calculation of entropy"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        e = calc.calculate_entropy()

        expected = {"gcf1": 0.0, "gcf2": 0.0}
        self.assertEqual(e, expected)

    def test_associations(self):
        """Test the calculation of associations"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        assoc = calc.compare_association()

        self.assertEqual(assoc, (1.0, 0.0, 1.0, 0.0))

    def test_conf_matrix(self):
        """Test the calculation confusion matrix"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        mat_data = calc.confusion_matrix()

        expected = ([[2, 0], [0, 1]], ["gcf1", "gcf2"], ["gcf1", "gcf2"])

        self.assertEqual(mat_data, expected)

    def test_summary(self):
        """Test the calculation of summary statistics"""
        curated = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}
        computed = {"bgc1": "gcf1", "bgc2": "gcf1", "bgc3": "gcf2"}

        calc = BenchmarkMetrics(curated, computed)
        summ = calc.calculate_summary()

        expected = (2, 2, 1.5, 1.5, 1, 1)
        self.assertEqual(summ, expected)
