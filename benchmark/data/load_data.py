"""Contains benchmarking functionality and evaluation"""

# from python
from pathlib import Path

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
        """Read tsv file with curated GCF assignments"""
        with open(self.curated_path) as inf:
            _ = inf.readline()
            data = {}
            for line in inf:
                bgc, family = line.strip().split("\t")
                data[bgc] = family
            self.curated_labels = data

    # TODO: loading runs with multiple GCF cutoffs
    def load_computed_labels(self):
        """Load computed GCF assignments from database"""
        DB.load_from_disk(self.db_path)
        fam_table = DB.metadata.tables["bgc_record_family"]
        record_table = DB.metadata.tables["bgc_record"]
        gbk_table = DB.metadata.tables["gbk"]

        select_query = (
            fam_table.join(record_table, fam_table.c.record_id == record_table.c.id)
            .join(gbk_table, record_table.c.gbk_id == gbk_table.c.id)
            .select()
            .add_columns(fam_table.c.record_id, fam_table.c.family, gbk_table.c.path)
            .compile()
        )

        cursor_results = DB.execute(select_query)
        self.computed_labels = {
            Path(result.path).stem: result.family for result in cursor_results
        }
