"""Contains benchmarking functionality and evaluation"""

# from python
from pathlib import Path
from collections import Counter

# from dependencies
from sklearn.metrics import v_measure_score
from scipy.stats import entropy

# from other modules
from big_scape.data import DB


class BenchmarkData:
    """Container for curated and computed GCF assignments

    Atrributes:
        curated_path (Path): location of curated GCF assignments file
        db_path (Path): location of BiG-SCAPE output database
    """

    def __init__(self, curated_path: Path, db_path: Path):
        self.curated_path = curated_path
        self.db_path = db_path

    def load_curated_labels(self):
        """Read tsv (?) file with curated GCF assignments"""
        with open(self.curated_path) as inf:
            _ = inf.readline()
            self.curated_labels = {
                bgc: family for line in inf for bgc, family in line.strip().split("\t")
            }

    def load_computed_labels(self):
        """Load copmuted GCF assignments from database"""
        DB.load_from_disk(self.db_path)
        fam_table = DB.metadata.tables["bgc_record_family"]

        select_query = (
            fam_table.select()
            .add_columns(fam_table.c.record_id, fam_table.c.family)
            .compile()
        )

        cursor_results = DB.execute(select_query)
        self.computed_labels = {
            result.record_id: result.family for result in cursor_results
        }

    def calculate_v_measure(self) -> float:
        """Calculate V-measure between curated and computed GCF assignments

        Following Rosenberg and Hirschberg:
        V-measure: A conditional entropy-based external cluster evaluation measure.
        """
        computed_bgcs = sorted(self.computed_labels.keys())

        curated_fams = [self.curated_labels[bgc] for bgc in computed_bgcs]
        computed_fams = [self.computed_labels[bgc] for bgc in computed_bgcs]

        return v_measure_score(curated_fams, computed_fams)

    def calculate_purity(self) -> float:
        """Calculate purity P of computed GCFs

        A score close to 1 indicates most computed clusters contains a single label based
        on the curated group of GCFs.
        """
        total_bgcs = len(self.computed_labels)

        computed_families: dict[str, list[str]] = {}
        for bgc, family in self.computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        max_label_occs: list[int] = []
        for family, bgcs in computed_families.items():
            curated_labels = [self.curated_labels[bgc] for bgc in bgcs]
            max_label_occs.append(Counter(curated_labels).most_common(1)[0][1])
        return sum(max_label_occs) / total_bgcs

    def calculate_entropy(self) -> dict[str, float]:
        """Calculate entropy for each computed GCF"""
        computed_families: dict[str, list[str]] = {}
        for bgc, family in self.computed_labels.items():
            computed_families.setdefault(family, []).append(bgc)

        entropies: dict[str, float] = {}
        for family, bgcs in computed_families.items():
            curated_labels = [self.curated_labels[bgc] for bgc in bgcs]
            pk = [curated_labels.count(label) for label in set(curated_labels)]
            entropies[family] = entropy(pk, base=2)
        return entropies
