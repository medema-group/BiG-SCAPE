"""Contains output generation and formatting"""

# from python
import os
import numpy as np
from pathlib import Path

# from dependencies
import matplotlib.pyplot as plt


class OutputGenerator:
    """Class to generate output files containing benchamrking results

    Attributes:
        output_dir (Path): Path pointing to base output directory
    """

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir

    def initialize_output_dir(self):
        """Set up output directory"""
        if not self.output_dir.exists():
            os.makedirs(self.output_dir)

    def output_purities(self, purities: dict[str, float]) -> None:
        """Write sorted computed GCF purities to output tsv file"""
        filename = self.output_dir / "GCF_purities.tsv"

        with open(filename, "w") as outf:
            outf.write("GCF name\tPurity\n")
            for family in sorted(purities, key=lambda f: purities[f], reverse=True):
                outf.write(f"{family}\t{purities[family]}\n")

    def output_entropies(self, entropies: dict[str, float]) -> None:
        """Write sorted computed GCF entropies to output tsv file"""
        filename = self.output_dir / "GCF_entropies.tsv"

        with open(filename, "w") as outf:
            outf.write("GCF name\tEntropy\n")
            for family in sorted(entropies, key=lambda f: entropies[f], reverse=True):
                outf.write(f"{family}\t{entropies[family]}\n")

    def output_summary(
        self,
        homogeneity: float,
        completeness: float,
        v_measure: float,
        purities: dict[str, float],
        entropies: dict[str, float],
        associations: tuple[float, float, float, float],
        summary_stats: tuple[float, float, float, float, int, int],
    ) -> None:
        """Output summary file comparing curated and computed GCF assignments

        Args:
            homogeneity: scores how homogeneous formed clusters are, ie a single label
            completeness: scores how complete curated families are in computed clusters
            v_measure: harmonic mean of homogeneity and completeness
            purities: similar to homogeneity, how pure are formed clusters
            entropies: describes the level of randomness in formed clusters
            associations: fraction of correct/wrong/present/missing assiciations per BGC
            summary_stats: number of families and singletons, average family size
        """
        filename = self.output_dir / "Summary.tsv"
        correct, wrong, present, missing = associations
        cur_fams, comp_fams, cur_size, comp_size, cur_sing, comp_sing = summary_stats
        with open(filename, "w") as outf:
            outf.write(
                "\tCurated_GCFs\tComputed_GCFs\n"
                + f"Number of families\t{cur_fams}\t{comp_fams}\n"
                + f"Number of singletons\t{cur_sing}\t{comp_sing}\n"
                + f"Average family size\t{cur_size:.2f}\t{comp_size:.2f}\n\n"
                + f"V-measure\t{v_measure:.4f}\n"
                + f"Homogeneity\t{homogeneity:.4f}\n"
                + f"Completeness\t{completeness:.4f}\n\n"
                + f"Average purity\t{np.average(list(purities.values())):.4f}\n"
                + f"Average entropy\t{np.average(list(entropies.values())):.4f}\n\n"
                + f"Fraction of correct associations per BGC\t{correct:.4f}\n"
                + f"Fraction of wrong   associations per BGC\t{wrong:.4f}\n"
                + f"Fraction of present curated associations per BGC\t{present:.4f}\n"
                + f"Fraction of missing curated associations per BGC\t{missing:.4f}\n"
            )

    def plot_per_cutoff(self, metrics):
        """Plot metrics per used cutoff"""
        cutoffs = metrics.keys()
        homogeneity = [metrics[cut]["homogeneity"] for cut in cutoffs]
        completeness = [metrics[cut]["completeness"] for cut in cutoffs]
        v_measure = [metrics[cut]["v_measure"] for cut in cutoffs]

        fig = plt.figure()
        ax = fig.gca()

        ax.plot(
            cutoffs,
            homogeneity,
            linewidth=0.5,
            linestyle="--",
            c="red",
            label="Homogeneity",
        )
        ax.plot(
            cutoffs,
            completeness,
            linewidth=0.5,
            linestyle="--",
            c="blue",
            label="Completeness",
        )
        ax.plot(cutoffs, v_measure, c="purple", label="V-measure")

        plt.title("External cluster evaluation metrics per used cutoff")
        plt.xlim(0, 1)
        plt.ylim(0, 1)
        plt.xlabel("BiG-SCAPE family cutoff")
        plt.ylabel("Score")
        plt.legend()
        plt.savefig(self.output_dir / "Scores_per_cutoff.png")
