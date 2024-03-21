"""Contains output generation and formatting"""

# from python
import os
import logging
import numpy as np
from pathlib import Path
from typing import Any, Optional

# from dependencies
import matplotlib.pyplot as plt

# suppress dependency debug logs
for logger in ["PIL", "matplotlib"]:
    logging.getLogger(logger).setLevel(logging.WARNING)


class OutputGenerator:
    """Class to generate output files containing benchamrking results

    Attributes:
        output_dir (Path): Path pointing to base output directory
        metadata (str): Metadata text to be writted in each output file
        run_name (str): Dataset name and starttime added to each filename
    """

    def __init__(self, output_dir: Path, metadata: str, name: str) -> None:
        self.output_dir = output_dir
        self.metadata = metadata
        self.name = name

    def initialize_output_dir(self) -> None:
        """Set up output directory"""
        if not self.output_dir.exists():
            os.makedirs(self.output_dir)

    @staticmethod
    def generate_metadata(run: dict, cutoff: Optional[str] = None) -> str:
        """Generate metadata based on the current comparison

        Args:
            run (dict): stored command line arguments
            cutoff (str): optional cutoff value being compared
        """
        t = run["start_time"]
        i = run["big_dir"]
        g = run["gcf_assignment_file"]
        meta = f"# {t}\n# Comparing input directory: {i}\n# To curated GCFs: {g}\n"
        if cutoff is not None:
            meta += f"# Used cutoff: {cutoff}\n"
        return meta + "\n"

    def output_metrics(
        self, metrics: dict[str, dict[str, Any]], fam_cutoff: str
    ) -> None:
        """Generates output for all calculated metrics"""
        self.initialize_output_dir()
        self.output_purities(metrics[fam_cutoff]["purities"])
        self.output_entropies(metrics[fam_cutoff]["entropies"])
        self.output_matrix(metrics[fam_cutoff]["conf_matrix"])
        self.output_summary(
            metrics[fam_cutoff]["homogeneity"],
            metrics[fam_cutoff]["completeness"],
            metrics[fam_cutoff]["v_measure"],
            metrics[fam_cutoff]["purities"],
            metrics[fam_cutoff]["entropies"],
            metrics[fam_cutoff]["associations"],
            metrics[fam_cutoff]["summary_stats"],
        )
        self.plot_conf_matrix_heatmap(metrics[fam_cutoff]["conf_matrix"])

    def output_purities(self, purities: dict[str, float]) -> None:
        """Write sorted computed GCF purities to output tsv file

        Args:
            purities: dict of purity per computed GCF
        """
        filename = self.output_dir / f"Purities_{self.name}.tsv"

        with open(filename, "w") as outf:
            outf.write(self.metadata)
            outf.write("GCF name\tPurity\n")
            for family in sorted(purities, key=lambda f: purities[f], reverse=True):
                outf.write(f"{family}\t{purities[family]}\n")

    def output_entropies(self, entropies: dict[str, float]) -> None:
        """Write sorted computed GCF entropies to output tsv file

        Args:
            entropies: dict of entropy per computed GCF
        """
        filename = self.output_dir / f"Entropies_{self.name}.tsv"

        with open(filename, "w") as outf:
            outf.write(self.metadata)
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
        filename = self.output_dir / f"Summary_{self.name}.tsv"
        correct, wrong, present, missing = associations
        cur_fams, comp_fams, cur_size, comp_size, cur_sing, comp_sing = summary_stats
        with open(filename, "w") as outf:
            outf.write(self.metadata)
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

    def output_matrix(
        self, matrix_data: tuple[list[list[int]], list[str], list[str]]
    ) -> None:
        """Write confusion matrix to tsv file

        Args:
            matrix_data: contains confusion matrix, row labels and column labels
        """
        filename = self.output_dir / f"Confusion_matrix_{self.name}.tsv"
        matrix, row_lab, col_lab = matrix_data
        with open(filename, "w") as outf:
            outf.write(self.metadata)
            col_lab_fmt = "\t".join(map(str, col_lab))
            outf.write(f"\t{col_lab_fmt}\n")
            for label, row in zip(row_lab, matrix):
                row_fmt = "\t".join(map(str, row))
                outf.write(f"{label}\t{row_fmt}\n")

    def plot_per_cutoff(self, metrics: dict[str, dict[str, Any]]) -> None:
        """Plot metrics per used cutoff

        Args:
            metrics: data dictionary storing all metrics per used cutoff
        """
        cutoffs = sorted(metrics.keys())
        homogeneity = [metrics[cut]["homogeneity"] for cut in cutoffs]
        completeness = [metrics[cut]["completeness"] for cut in cutoffs]
        v_measure = [metrics[cut]["v_measure"] for cut in cutoffs]
        wrong = [metrics[cut]["associations"][1] for cut in cutoffs]
        missing = [metrics[cut]["associations"][3] for cut in cutoffs]
        cutoffs_fl = list(map(float, cutoffs))

        fig = plt.figure()
        ax = fig.gca()  # type: ignore

        h = ax.plot(
            cutoffs_fl,
            homogeneity,
            linestyle="--",
            c="#FFC107",  # yellow
            label="Homogeneity",
        )
        c = ax.plot(
            cutoffs_fl,
            completeness,
            linestyle="--",
            c="#1E88E5",  # blue
            label="Completeness",
        )
        v = ax.plot(cutoffs_fl, v_measure, c="black", label="V-measure")

        wl = ax.plot(
            cutoffs_fl,
            wrong,
            linestyle="-.",
            c="#D81B60",  # red
            label="Wrong links",
        )
        ml = ax.plot(
            cutoffs_fl,
            missing,
            linestyle="-.",
            c="#E68981",  # pink
            label="Missing links",
        )
        ax.text(0, -0.4, self.metadata, transform=ax.transAxes)
        plt.title("External cluster evaluation metrics per used cutoff")
        plt.xlabel("BiG-SCAPE family cutoff")
        plt.ylabel("Score")
        plots = h + c + v + wl + ml
        ax.legend(plots, [p.get_label() for p in plots], loc=0)
        plt.savefig(
            self.output_dir / f"Scores_per_cutoff_{self.name}.png",
            bbox_inches="tight",
            dpi=400,
        )

    def plot_conf_matrix_heatmap(
        self, matrix_data: tuple[list[list[int]], list[str], list[str]]
    ) -> None:
        """Plot confusion matrix as heatmap

        Args:
            matrix_data: contains confusion matrix, row labels and column labels
        """
        matrix, row_lab, col_lab = matrix_data
        plt.imshow(matrix, cmap="binary", interpolation=None)  # type: ignore
        ax = plt.gca()  # type: ignore
        ax.set_xticks(
            range(len(col_lab)),
            labels=col_lab,
            fontsize=1,
            rotation=90,
            rotation_mode="default",
            ha="center",
        )
        ax.set_yticks(range(len(row_lab)), labels=row_lab, fontsize=1)
        plt.xlabel("Computed GCFs")
        plt.ylabel("Curated GCFs")
        plt.title("Overlap of curated and computed GCFs")
        plt.savefig(
            self.output_dir / f"Confusion_heatmap_{self.name}.png",
            bbox_inches="tight",
            dpi=700,
        )

    def output_summary_per_cutoff(self, metrics: dict[str, dict[str, Any]]) -> None:
        """Write overview of evaluation metrics across all used cutoffs

        Args:
            metrics: data dictionary storing all metrics per used cutoff
        """
        filename = self.output_dir / f"Benchmark_summary_{self.name}.txt"

        with open(filename, "w") as outf:
            outf.write(self.metadata)
            cutoff_fmt = "\t".join(list(metrics.keys()))

            v_fmt = "\t".join([f"{metrics[c]['v_measure']:.3f}" for c in metrics])
            h_fmt = "\t".join([f"{metrics[c]['homogeneity']:.3f}" for c in metrics])
            c_fmt = "\t".join([f"{metrics[c]['completeness']:.3f}" for c in metrics])

            purity = [np.mean(list(metrics[c]["purities"].values())) for c in metrics]
            p_fmt = "\t".join([f"{p:.3f}" for p in purity])
            entropy = [np.mean(list(metrics[c]["entropies"].values())) for c in metrics]
            e_fmt = "\t".join(f"{e:.3f}" for e in entropy)

            cur_fam_nr = [metrics[c]["summary_stats"][0] for c in metrics]
            comp_fam_nr = [metrics[c]["summary_stats"][1] for c in metrics]
            fam_nr_diffs = [comp - cur for cur, comp in zip(cur_fam_nr, comp_fam_nr)]
            nr_diff_fmt = "\t".join(map(str, fam_nr_diffs))

            outf.write(
                f"Used_cutoff\t{cutoff_fmt}\n"
                + f"V-measure\t{v_fmt}\n"
                + f"Homogeneity\t{h_fmt}\n"
                + f"Completeness\t{c_fmt}\n"
                + f"Purity\t{p_fmt}\n"
                + f"Entropy\t{e_fmt}\n"
                + f"GCF number diff\t{nr_diff_fmt}\n"
            )
