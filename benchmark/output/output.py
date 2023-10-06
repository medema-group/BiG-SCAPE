"""Contains output generation and formatting"""

# from python
from pathlib import Path


class OutputGenerator:
    """Class to generate output files containing benchamrking results"""

    def __init__(self, output_dir: Path) -> None:
        self.output_dir = output_dir

    def output_entropies(self, entropies: dict[str, float]) -> None:
        """Write sorted computed GCF entropies to output tsv file"""
        filename = self.output_dir / "GCF_entropies.tsv"

        with open(filename, "w") as outf:
            outf.write("GCF name\tEntropy\n")
            for family in sorted(entropies, key=lambda f: entropies[f], reverse=True):
                outf.write(f"{family}\t{entropies[family]}\n")
