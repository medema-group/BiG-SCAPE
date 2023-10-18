"""Contains benchmark data loading and storage"""

# from python
from pathlib import Path
from datetime import datetime
import sqlite3

# from other modules
from big_scape.data import DB


class BenchmarkData:
    """Container for curated and computed GCF assignments

    Atrributes:
        curated_path (Path): location of curated GCF assignments file
        in_path (Path): location of BiG-SCAPE output database
    """

    def __init__(self, curated_path: Path, in_path: Path) -> None:
        self.curated_path = curated_path
        self.in_path = in_path

    def read_gcf_tsv(self, filename: Path) -> dict[str, str]:
        """Read GCF assignments tsv file

        Expects tsv file with BGC name and family per line and stores in {bgc_name: fam}

        Args:
            filename (Path): tsv file to be read

        Returns:
            data dictionary with BGC name linking to family assignment
        """
        with open(filename) as inf:
            _ = inf.readline()
            data = {}
            for line in inf:
                bgc, family = line.strip().split("\t")
                data[bgc] = family
            return data

    def load_curated_labels(self) -> None:
        """Read tsv file with curated GCF assignments"""
        self.curated_labels = self.read_gcf_tsv(self.curated_path)

    def load_computed_labels(self) -> None:
        """Load computed GCF assignments from BS1/BS2/BSLICE output

        Raises:
            FileNotFoundError: Unexpected structure of BS1/BS2/BSLICE output directories
        """
        bs2_path = self.in_path / "data_sqlite.db"
        bs1_path = self.in_path / "network_files"
        bslice_path = self.in_path / "result/data.db"
        if bs2_path.exists():
            self.load_computed_bs2_labels(bs2_path)
        elif bs1_path.exists():
            self.load_computed_bs1_labels(bs1_path)
        elif bslice_path.exists():
            self.load_computed_bslice_labels(bslice_path)
        else:
            raise FileNotFoundError("Unable to locate correct output files")

    def load_computed_bs2_labels(self, db_path: Path) -> None:
        """Load computed GCF assignments from BS2 database to {cutoff: {bgc_name: fam}}

        Args:
            db_path (Path): Path to existing BS2 database
        """
        DB.load_from_disk(db_path)
        if DB.metadata is not None:
            fam_table = DB.metadata.tables["bgc_record_family"]
            record_table = DB.metadata.tables["bgc_record"]
            gbk_table = DB.metadata.tables["gbk"]
        else:
            raise RuntimeError("DB.metadata is None")

        select_query = (
            fam_table.join(record_table, fam_table.c.record_id == record_table.c.id)
            .join(gbk_table, record_table.c.gbk_id == gbk_table.c.id)
            .select()
            .add_columns(fam_table.c.cutoff, fam_table.c.family, gbk_table.c.path)
            .compile()
        )

        cursor_results = DB.execute(select_query)
        self.computed_labels: dict[float, dict[str, str]] = {}
        for result in cursor_results:
            bgc_name = str(Path(result.path).stem)
            family = str(result.family)
            self.computed_labels.setdefault(result.cutoff, {})[bgc_name] = family

        # add missing singletons per cutoff, assign record_id as family id
        select_query = (
            record_table.join(gbk_table, record_table.c.gbk_id == gbk_table.c.id)
            .select()
            .add_columns(gbk_table.c.path, record_table.c.id)
            .compile()
        )
        cursor_results = DB.execute(select_query)
        name_to_record_id = {
            str(Path(result.path).stem): str(result.id) for result in cursor_results
        }

        for cutoff in self.computed_labels.keys():
            for name in name_to_record_id.keys():
                if name not in self.computed_labels[cutoff]:
                    self.computed_labels[cutoff][name] = name_to_record_id[name]

    def load_computed_bs1_labels(self, data_path: Path) -> None:
        """Load computed GCFs from BS1 results files to {cutoff: {bgc_name: fam}}

        Args:
            data_path (Path): Path pointing to network files of BS1 output directory

        Raises:
            FileNotFoundError: Missing BS1 results in given output directory
        """
        runs = list(data_path.glob("*"))
        if len(runs) == 0:
            raise FileNotFoundError("No BiG-SCAPE 1 output found")
        elif len(runs) == 1:
            run_path = runs.pop(0) / "mix"
        else:
            # select the most recent run
            most_recent = datetime(1, 1, 1)
            for run in runs:
                date, time, _ = str(run.stem).split("_", 2)
                year, month, day = map(int, date.split("-"))
                hour, minute, second = map(int, time.split("-"))
                run_datetime = datetime(year, month, day, hour, minute, second)
                if run_datetime > most_recent:
                    most_recent = run_datetime
                    run_path = run / "mix"

        if not run_path.exists():
            raise FileNotFoundError("No BiG-SCAPE mix results found in most recent run")

        self.computed_labels = {}
        for clustering_file in run_path.glob("*.tsv"):
            cutoff = float(clustering_file.stem.rpartition("c")[-1])
            self.computed_labels[cutoff] = self.read_gcf_tsv(clustering_file)

    def load_computed_bslice_labels(self, db_path: Path) -> None:
        """Load computed GCFs from BiG-SLICE output to {threshold: {bgc: family}}

        Args:
            db_path (Path): Path pointing to BiG-SLICE output database
        """
        con = sqlite3.connect(db_path)
        cur = con.cursor()

        # if multiple runs are present with the same threshold, only take most recent
        thresh_data = cur.execute(
            "SELECT clustering.threshold, clustering.id FROM clustering"
        )
        threshs = {thresh: run_id for thresh, run_id in thresh_data}

        # collect bgc and their family assignment per threshold
        cursor_results = cur.execute(
            "SELECT bgc.orig_filename, gcf_membership.gcf_id, clustering.threshold "
            "FROM bgc "
            "INNER JOIN gcf_membership ON bgc.id==gcf_membership.bgc_id "
            "INNER JOIN gcf ON gcf.id==gcf_membership.gcf_id "
            "INNER JOIN clustering ON clustering.id==gcf.clustering_id "
            f"WHERE clustering.id IN ({','.join(map(str, threshs.values()))})"
        )
        self.computed_labels = {}
        for result in cursor_results:
            bgc_name = str(Path(result[0]).stem)
            family = str(result[1])
            thresh = result[2]
            self.computed_labels.setdefault(thresh, {})[bgc_name] = family
