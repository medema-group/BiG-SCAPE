"""Contains benchmarking functionality and evaluation"""

# from python
from pathlib import Path

# from dependencies
from sklearn.metrics import v_measure_score

# from other modules
from big_scape.data import DB


class BenchmarkData:
    """Container for curated and computed GCF assignments

    Atrributes:

    """

    def __init__(self, curated_gcfs: Path, db_path: Path):
        self.curated_gcfs = curated_gcfs
        self.db_path = db_path

    def load_curated_labels(self):
        """Read tsv (?) file with curated GCF assignments"""
        with open(self.curated_gcfs) as inf:
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

    def calculate_v_measure(self):
        """Calculate V-measure between curated and computed GCF assignments

        Following Rosenberg and Hirschberg:
        V-measure: A conditional entropy-based external cluster evaluation measure.
        """
        computed_bgcs = sorted(self.computed_labels.keys())

        curated_fams = [self.curated_labels[bgc] for bgc in computed_bgcs]
        computed_fams = [self.computed_labels[bgc] for bgc in computed_bgcs]

        return v_measure_score(curated_fams, computed_fams)
