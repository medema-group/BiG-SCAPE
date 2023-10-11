"""Contains benchmark data loading and storage"""

# from python
from pathlib import Path
from datetime import datetime

# from other modules
from big_scape.data import DB


class BenchmarkData:
    """Container for curated and computed GCF assignments

    Atrributes:
        curated_path (Path): location of curated GCF assignments file
        db_path (Path): location of BiG-SCAPE output database
    """

    def __init__(self, curated_path: Path, in_path: Path):
        self.curated_path = curated_path
        self.in_path = in_path

    def read_gcf_tsv(self, filename):
        """Read GCF assignments tsv file

        Expects tsv file with BGC name and family per line and stores in {bgc_name: fam}

        Args:
            filename (Path): tsv file to be read
        """
        with open(filename) as inf:
            _ = inf.readline()
            data = {}
            for line in inf:
                bgc, family = line.strip().split("\t")
                data[bgc] = family
            return data

    def load_curated_labels(self):
        """Read tsv file with curated GCF assignments"""
        self.curated_labels = self.read_gcf_tsv(self.curated_path)

    def load_computed_labels(self):
        """Load computed GCF assignments from BS1/BS2/BSLICE output"""
        bs2_path = self.in_path / "data_sqlite.db"
        bs1_path = self.in_path / "network_files"
        bsl_path = self.in_path / "TBD"
        if bs2_path.exists():
            self.load_computed_bs2_labels(bs2_path)
        elif bs1_path.exists():
            self.load_computed_bs1_labels(bs1_path)
        elif bsl_path.exists():
            self.load_computed_bsl_labels(bsl_path)
        else:
            raise FileNotFoundError("Unable to locate correct ouput files")

    def load_computed_bs2_labels(self, db_path):
        """Load computed GCF assignments from BS2 database to {cutoff: {bgc_name: fam}}

        Args:
            db_path (Path): Path to existing BS2 database
        """
        DB.load_from_disk(db_path)
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

    def load_computed_bs1_labels(self, data_path):
        """Load computed GCFs from BS1 results files to {cutoff: {bgc_name: fam}}

        Args:
            data_path (Path): Path pointing to network files of BS1 output directory
        """
        runs = list(data_path.glob("*"))
        if len(runs) == 0:
            raise FileNotFoundError("No BiG-SCAPE 1 output found")
        elif len(runs) == 1:
            run_path = data_path / runs.pop(0) / "mix"
        else:
            # select the most recent run
            most_recent = datetime(1, 1, 1)
            for run in runs:
                date, time, _ = str(run.stem).split("_", 2)
                run_datetime = datetime(
                    *map(int, date.split("-")), *map(int, time.split("-"))
                )
                if run_datetime > most_recent:
                    most_recent = run_datetime
                    run_path = run / "mix"

        if not run_path.exists():
            raise FileNotFoundError("No BiG-SCAPE mix results found in most recent run")

        self.computed_labels = {}
        for clustering_file in run_path.glob("*.tsv"):
            cutoff = float(clustering_file.stem.rpartition("c")[-1])
            self.computed_labels[cutoff] = self.read_gcf_tsv(clustering_file)

    def load_computed_bsl_labels(self, data_path):
        """Load computed GCFs from BiG-SLICE output to {threshold: {bgc: family}}

        Args:
            data_path (Path): Path pointing to BiG-SLICE output directory
        """
        # TODO: parse BiG-SLICE output
        raise NotImplementedError()
