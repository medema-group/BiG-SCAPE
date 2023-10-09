"""Contains output generation and formatting"""

# from python
import os
from pathlib import Path


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

    def output_v_measure(self, v_measure: float) -> None:
        """Write computed GCF V measure to output txt file"""
        filename = self.output_dir / "GCF_V_measure.txt"

        with open(filename, "w") as outf:
            outf.write(f"Computed V-measure: {v_measure}\n")

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
