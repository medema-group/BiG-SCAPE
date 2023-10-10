"""Contains benchmark data loading and storage"""

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
        """Read tsv file with curated GCF assignments

        Expects tsv file with BGC name and family per line and stores in {bgc_name: fam}
        """
        with open(self.curated_path) as inf:
            _ = inf.readline()
            data = {}
            for line in inf:
                bgc, family = line.strip().split("\t")
                data[bgc] = family
            self.curated_labels = data

    def load_computed_labels(self):
        """Load computed GCF assignments from database

        Parse an existing BiG-SCAPE database into {cutoff: {bgc_name: fam}}
        """
        DB.load_from_disk(self.db_path)
        fam_table = DB.metadata.tables["bgc_record_family"]
        record_table = DB.metadata.tables["bgc_record"]
        gbk_table = DB.metadata.tables["gbk"]

        select_query = (
            fam_table.join(record_table, fam_table.c.record_id == record_table.c.id)
            .join(gbk_table, record_table.c.gbk_id == gbk_table.c.id)
            .select()
            .add_columns(fam_table.c.cutoff, fam_table.c.family, gbk_table.c.path)
            .compile()
        )

        cursor_results = DB.execute(select_query)
        self.computed_labels = {}
        for result in cursor_results:
            bgc_name = Path(result.path).stem
            self.computed_labels.setdefault(result.cutoff, {})[bgc_name] = result.family

        # add missing singletons per cutoff, assign record_id as family id
        select_query = (
            record_table.join(gbk_table, record_table.c.gbk_id == gbk_table.c.id)
            .select()
            .add_columns(gbk_table.c.path, record_table.c.id)
            .compile()
        )
        cursor_results = DB.execute(select_query)
        name_to_record_id = {
            Path(result.path).stem: result.id for result in cursor_results
        }

        for cutoff in self.computed_labels.keys():
            for name in name_to_record_id.keys():
                if name not in self.computed_labels[cutoff]:
                    self.computed_labels[cutoff][name] = name_to_record_id[name]
