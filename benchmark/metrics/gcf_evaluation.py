"""Contains calculation of benchmarking evaluation metrics"""

# from python
from collections import Counter

# from dependencies
from sklearn.metrics import v_measure_score
from scipy.stats import entropy


class BenchmarkMetrics:
    """Class to store evaluation metric calculation functions"""

    @staticmethod
    def calculate_v_measure(
        curated_labels: dict[str, str], computed_labels: dict[str, str]
    ) -> float:
        """Calculate V-measure between curated and computed GCF assignments

        Following Rosenberg and Hirschberg:
        V-measure: A conditional entropy-based external cluster evaluation measure.

        Args:
            curated_labels (dict[str, str]): curated GCF assignments per BGC
            computed_labels (dict[str, str]): BiG-SCAPE computed GCF assignments per BGC
        Returns:
            float: V-measure
        """
        computed_bgcs = computed_labels.keys()

        curated_fams = [curated_labels[bgc] for bgc in computed_bgcs]
        computed_fams = [computed_labels[bgc] for bgc in computed_bgcs]

        return v_measure_score(curated_fams, computed_fams)

    @staticmethod
    def calculate_purity(
        curated_labels: dict[str, str], computed_labels: dict[str, str]
    ) -> dict[str, float]:
        """Calculate purity P of each computed GCF

        A score close to 1 indicates most computed clusters contains a single label based
        on the curated group of GCFs.

        Args:
            curated_labels (dict[str, str]): curated GCF assignments per BGC
            computed_labels (dict[str, str]): BiG-SCAPE computed GCF assignments per BGC
        Returns:
            dict[str, float]: dictionary containing each family and their purity
        """
        computed_families: dict[str, list[str]] = {}
        for bgc, family in computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        purities: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels_in_fam = [curated_labels[bgc] for bgc in bgcs]
            max_label_occ = Counter(curated_labels_in_fam).most_common(1)[0][1]
            purities[family] = max_label_occ / len(curated_labels_in_fam)
        return purities

    @staticmethod
    def calculate_entropy(
        curated_labels: dict[str, str], computed_labels: dict[str, str]
    ) -> dict[str, float]:
        """Calculate entropy of each computed GCF

        Args:
            curated_labels (dict[str, str]): curated GCF assignments per BGC
            computed_labels (dict[str, str]): BiG-SCAPE computed GCF assignments per BGC
        Returns:
            dict[str, float]: dictionary containing each family and their entropy
        """
        computed_families: dict[str, list[str]] = {}
        for bgc, family in computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        entropies: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels_in_fam = [curated_labels[bgc] for bgc in bgcs]
            pk = [curated_labels_in_fam.count(label) for label in set(curated_labels)]
            entropies[family] = entropy(pk, base=2)
        return entropies
