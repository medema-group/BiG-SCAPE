"""Contains calculation of benchmarking evalutaion metrics"""

# from python
from collections import Counter

# from dependencies
from sklearn.metrics import v_measure_score
from scipy.stats import entropy

# from other modules
from benchmark.data import BenchmarkData


class BenchmarkMetrics:
    """Class to store evaluation metrics calculation functions"""

    @staticmethod
    def calculate_v_measure(data: BenchmarkData) -> float:
        """Calculate V-measure between curated and computed GCF assignments

        Following Rosenberg and Hirschberg:
        V-measure: A conditional entropy-based external cluster evaluation measure.
        """
        computed_bgcs = sorted(data.computed_labels.keys())

        curated_fams = [data.curated_labels[bgc] for bgc in computed_bgcs]
        computed_fams = [data.computed_labels[bgc] for bgc in computed_bgcs]

        return v_measure_score(curated_fams, computed_fams)

    @staticmethod
    def calculate_purity(data: BenchmarkData) -> float:
        """Calculate purity P of computed GCFs

        A score close to 1 indicates most computed clusters contains a single label based
        on the curated group of GCFs.
        """
        total_bgcs = len(data.computed_labels)

        computed_families: dict[str, list[str]] = {}
        for bgc, family in data.computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        max_label_occs: list[int] = []
        for family, bgcs in computed_families.items():
            curated_labels = [data.curated_labels[bgc] for bgc in bgcs]
            max_label_occs.append(Counter(curated_labels).most_common(1)[0][1])
        return sum(max_label_occs) / total_bgcs

    @staticmethod
    def calculate_entropy(data: BenchmarkData) -> dict[str, float]:
        """Calculate entropy for each computed GCF"""
        computed_families: dict[str, list[str]] = {}
        for bgc, family in data.computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        entropies: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels = [data.curated_labels[bgc] for bgc in bgcs]
            pk = [curated_labels.count(label) for label in set(curated_labels)]
            entropies[family] = entropy(pk, base=2)
        return entropies
