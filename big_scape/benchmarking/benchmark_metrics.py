"""Contains calculation of benchmarking evaluation metrics"""

# from python
from collections import Counter
from typing import Any

# from dependencies
import numpy as np
from sklearn.metrics import homogeneity_completeness_v_measure
from scipy.stats import entropy


class BenchmarkMetrics:
    """Class to store evaluation metric calculation functions

    Attributes:
        curated_labels (dict[str, str]): curated GCF assignments for each BGC
        computed_labels (dict[str, str]): computed GCF assignments for each BGC
    """

    def __init__(self, curated_labels: dict[str, str], computed_labels: dict[str, str]):
        self.curated_labels = curated_labels
        self.computed_labels = computed_labels

    def calculate_metrics(self) -> dict[str, Any]:
        """Calculate all metrics

        Includes: homogeneity, completeness, v_measure, purities, entropies,
        conf_matrix, associations and summary_stats

        Returns:
            Dictionary holding all calculated metrics
        """
        homogeneity, completeness, v_measure = self.calculate_v_measure()
        purities = self.calculate_purity()
        entropies = self.calculate_entropy()
        conf_matrix = self.confusion_matrix()
        associations = self.compare_association()
        summary_stats = self.calculate_summary()
        return {
            "homogeneity": homogeneity,
            "completeness": completeness,
            "v_measure": v_measure,
            "purities": purities,
            "entropies": entropies,
            "conf_matrix": conf_matrix,
            "associations": associations,
            "summary_stats": summary_stats,
        }

    @staticmethod
    def create_fam_index(labels: dict[str, str]) -> dict[str, list[str]]:
        """Create family index linking each family to their BGC members

        Args:
            labels (dict[str, str]): dict linking each BGC to their assigned family

        Returns:
            Dictionary linking each family to a list of their BGC members
        """
        fam_index: dict[str, list[str]] = {}
        for bgc, fam in labels.items():
            fam_index.setdefault(fam, []).append(bgc)
        return fam_index

    def calculate_v_measure(self) -> tuple[float, float, float]:
        """Calculate V-measure between curated and computed GCF assignments

        Following Rosenberg and Hirschberg:
        V-measure: A conditional entropy-based external cluster evaluation measure.

        Returns:
            tuple[float, float, float]: homogeneity, completeness and V-measure
        """
        computed_bgcs = self.computed_labels.keys()

        curated_fams = [self.curated_labels[bgc] for bgc in computed_bgcs]
        computed_fams = [self.computed_labels[bgc] for bgc in computed_bgcs]

        return homogeneity_completeness_v_measure(curated_fams, computed_fams)

    def calculate_purity(self) -> dict[str, float]:
        """Calculate purity P of each computed GCF

        A score close to 1 indicates most computed clusters contains a single label based
        on the curated group of GCFs.

        Returns:
            dict[str, float]: dictionary containing each family and their purity
        """
        computed_families = BenchmarkMetrics.create_fam_index(self.computed_labels)

        purities: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels_in_fam = [self.curated_labels[bgc] for bgc in bgcs]
            max_label_occ = Counter(curated_labels_in_fam).most_common(1)[0][1]
            purities[family] = max_label_occ / len(curated_labels_in_fam)
        return purities

    def calculate_entropy(self) -> dict[str, float]:
        """Calculate entropy of each computed GCF

        Returns:
            dict[str, float]: dictionary containing each family and their entropy
        """
        computed_families = BenchmarkMetrics.create_fam_index(self.computed_labels)

        entropies: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels_in_fam = [self.curated_labels[bgc] for bgc in bgcs]
            pk = [
                curated_labels_in_fam.count(label)
                for label in set(curated_labels_in_fam)
            ]
            entropies[family] = entropy(pk, base=2)
        return entropies

    def compare_association(self) -> tuple[float, float, float, float]:
        """Collect fractions of correct/wrong and present/missing associations per BGC

        Returns:
            Fractions of correct and wrong associations per BGC based on computed
            families, and fraction of present and missing associations based on curated
            families.
        """
        computed_bgcs = self.computed_labels.keys()
        curated_labels_used = {
            bgc: fam for bgc, fam in self.curated_labels.items() if bgc in computed_bgcs
        }
        computed_fams = BenchmarkMetrics.create_fam_index(self.computed_labels)
        curated_fams = BenchmarkMetrics.create_fam_index(curated_labels_used)

        correct = []
        wrong = []
        present = []
        missing = []
        for bgc in computed_bgcs:
            cur_assoc = set(curated_fams[curated_labels_used[bgc]])
            comp_assoc = set(computed_fams[self.computed_labels[bgc]])
            # remove self-association
            cur_assoc.remove(bgc)
            comp_assoc.remove(bgc)

            # in the case of singletons, no associations are found
            if comp_assoc:
                correct.append(
                    len(cur_assoc.intersection(comp_assoc)) / len(comp_assoc)
                )
                wrong.append(len(comp_assoc - cur_assoc) / len(comp_assoc))
            if cur_assoc:
                present.append(len(cur_assoc.intersection(comp_assoc)) / len(cur_assoc))
                missing.append(
                    len([m for m in cur_assoc if m not in comp_assoc]) / len(cur_assoc)
                )

        return (
            float(np.average(correct)),
            float(np.average(wrong)),
            float(np.average(present)),
            float(np.average(missing)),
        )

    def confusion_matrix(self) -> tuple[list[list[int]], list[str], list[str]]:
        """Confusion matrix showing member overlap between curated and computed families

        Rows represent curated GCF assignments, columns computed GCF assignments. Matrix
        is sorted based on highest overlap between curated and computed families

        Returns:
            A sorted confusion matrix and its row and column GCF labels
        """
        mat: list[list[int]] = []

        comp_bgcs = self.computed_labels.keys()
        curated_labels_used = {
            bgc: fam for bgc, fam in self.curated_labels.items() if bgc in comp_bgcs
        }
        computed_fams = BenchmarkMetrics.create_fam_index(self.computed_labels)
        curated_fams = BenchmarkMetrics.create_fam_index(curated_labels_used)

        row_labels = sorted(curated_fams.keys())
        col_labels = sorted(computed_fams.keys())
        for cur_fam in row_labels:
            row: list[int] = []
            cur_members = set(curated_fams[cur_fam])
            for comp_fam in col_labels:
                comp_members = set(computed_fams[comp_fam])
                row.append(len(cur_members.intersection(comp_members)))
            mat.append(row)

        sorted_mat: list[list[int]] = []
        sorted_row_lab: list[str] = []
        for idx in range(len(col_labels)):
            col_vals = np.array([row[idx] for row in mat])
            max_rows = np.flatnonzero(col_vals == max(col_vals))
            for max_idx in max_rows:
                if row_labels[max_idx] not in sorted_row_lab:
                    sorted_mat.append(mat[max_idx])
                    sorted_row_lab.append(row_labels[max_idx])
        return sorted_mat, sorted_row_lab, col_labels

    def calculate_summary(self) -> tuple[int, int, float, float, int, int]:
        """Calculate summary GCF number and size for curated and computed GCFs

        Returns:
            Tuple containing number of curated and computed GCFs and their average size
        """
        computed_bgcs = self.computed_labels.keys()

        curated_labels_used = {
            bgc: fam for bgc, fam in self.curated_labels.items() if bgc in computed_bgcs
        }

        curated_fams = set(curated_labels_used.values())
        computed_fams = set(self.computed_labels.values())

        curated_avg_size = float(
            np.average(
                [list(curated_labels_used.values()).count(fam) for fam in curated_fams]
            )
        )
        computed_avg_size = float(
            np.average(
                [
                    list(self.computed_labels.values()).count(fam)
                    for fam in computed_fams
                ]
            )
        )
        curated_single = list(Counter(curated_labels_used.values()).values()).count(1)
        computed_single = list(Counter(self.computed_labels.values()).values()).count(1)
        return (
            len(curated_fams),
            len(computed_fams),
            curated_avg_size,
            computed_avg_size,
            curated_single,
            computed_single,
        )
