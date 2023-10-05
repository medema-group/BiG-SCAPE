"""Contains benchmarking functionality and evaluation"""

# from python
from pathlib import Path

# from dependencies

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
        computed_labels = {}
        for result in cursor_results:
            if result.family not in computed_labels:
                computed_labels[result.family] = set([result.record_id])
            else:
                computed_labels[result.family].add(result.record_id)
        self.computed_labels = computed_labels
