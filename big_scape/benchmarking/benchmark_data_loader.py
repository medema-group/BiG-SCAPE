"""Contains benchmark data loading and storage"""

# from python
from pathlib import Path
from datetime import datetime
from typing import TextIO
import sqlite3
import logging


class BenchmarkData:
    """Container for curated and computed GCF assignments

    Atrributes:
        curated_path (Path): location of curated GCF assignments file
        in_path (Path): location of BiG-SCAPE output database
    """

    def __init__(self, curated_path: Path, in_path: Path) -> None:
        self.curated_path = curated_path
        self.in_path = in_path

    def read_gcf_short_tsv(self, infile: TextIO) -> dict[str, str]:
        """Read GCF assignments tsv file short format

        Expects tsv file with BGC name and family per line and stores in {bgc_name: fam}

        Args:
            infile (TextIO): opened tsv file object

        Returns:
            dict[str, str]: data dictionary with BGC name linking to family assignment
        """
        data = {}
        for line in infile:
            bgc, family = line.strip().split("\t")
            data[bgc] = family
        return data

    def read_gcf_long_tsv(self, infile: TextIO) -> dict[str, str]:
        """Read GCF assignment tsv file long format

        Expects tsv file with BGC name, record type, record number and family per line.
        Stores record information as BGC_record_number, if record type is region ignores
        region number.

        Args:
            filename (TextIO): opened tsv file object

        Returns:
            dict[str, str]: dictionary linking BGC record to family
        """
        data = {}
        for line in infile:
            parts = line.strip().split("\t")
            clean_name = Path(parts[0]).name.replace(".gbk", "")

            if parts[1] == "region":
                bgc = clean_name
            else:
                bgc = f"{clean_name}_{parts[1]}_{parts[2]}"

            data[bgc] = parts[3]
        return data

    def load_curated_labels(self) -> None:
        """Read tsv file with curated GCF assignments"""
        logging.debug("Loading curated GCFs: %s", self.curated_path)

        with open(self.curated_path) as inf:
            header = inf.readline()

            if len(header.strip().split("\t")) == 2:
                logging.debug("Loading as short format tsv")
                self.curated_labels = self.read_gcf_short_tsv(inf)
            else:
                logging.debug("Loading as long format tsv")
                self.curated_labels = self.read_gcf_long_tsv(inf)

        if len(self.curated_labels) == 0:
            logging.warning("GCF assignment file is empty: %s", self.curated_path)
        logging.debug("Loaded %s BGC <-> GCF assignments", len(self.curated_labels))

    def load_computed_labels(self) -> None:
        """Load computed GCF assignments from BS1/BS2/BSLICE output

        Raises:
            FileNotFoundError: Unexpected structure of BS1/BS2/BSLICE output directories
        """
        bs2_path = self.in_path / "output_files"
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

    def load_computed_bs2_labels(self, data_path: Path) -> None:
        """Load computed GCFs from BS2 output files to {cutoff: {bgc_name: family}}

        Args:
            data_path (Path): Path pointing to output files of BS2 output directory
        """
        logging.info("Loading computed GCFs from BiG-SCAPE 2 output")
        run_times = [
            p.stem.replace("_full", "") for p in data_path.glob("*_full.network")
        ]
        if len(run_times) == 0:
            raise FileNotFoundError("No BiG-SCAPE 2 output found")
        elif len(run_times) == 1:
            run_datetime = run_times.pop(0)
        else:
            # select most recent run
            run_time = datetime(1, 1, 1)
            for dt_str in run_times:
                date, time = dt_str.split("_")
                day, month, year = map(int, date.split("-"))
                hour, minute, second = map(int, time.split("-"))
                current_dt = datetime(year, month, day, hour, minute, second)
                if current_dt > run_time:
                    run_time = current_dt
                    run_datetime = dt_str

        self.computed_labels = {}
        for clustering_file in data_path.glob(f"{run_datetime}*/mix/*clustering_*.tsv"):
            cutoff = clustering_file.stem.rpartition("_c")[-1]
            with open(clustering_file) as inclust:
                inclust.readline()
                self.computed_labels[cutoff] = self.read_gcf_long_tsv(inclust)

        logging.info(
            "Found clusterings containing %s BGCs across %s cutoffs",
            len(self.computed_labels[cutoff]),
            len(self.computed_labels),
        )

    def load_computed_bs1_labels(self, data_path: Path) -> None:
        """Load computed GCFs from BS1 results files to {cutoff: {bgc_name: fam}}

        Args:
            data_path (Path): Path pointing to network files of BS1 output directory

        Raises:
            FileNotFoundError: Missing BS1 results in given output directory
        """
        logging.info("Loading computed GCFs from BiG-SCAPE 1 output")
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
        for clustering_file in run_path.glob("*_clustering_*.tsv"):
            cutoff = clustering_file.stem.rpartition("_c")[-1]
            with open(clustering_file) as inclust:
                inclust.readline()
                self.computed_labels[cutoff] = self.read_gcf_short_tsv(inclust)

        logging.info(
            "Found clusterings containing %s BGCs across %s cutoffs",
            len(self.computed_labels[cutoff]),
            len(self.computed_labels),
        )

    def load_computed_bslice_labels(self, db_path: Path) -> None:
        """Load computed GCFs from BiG-SLICE output to {threshold: {bgc: family}}

        Args:
            db_path (Path): Path pointing to BiG-SLICE output database
        """
        logging.info("Loading computed GCFs from BiG-SLICE output")
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
            thresh = str(result[2])
            self.computed_labels.setdefault(thresh, {})[bgc_name] = family

        logging.info(
            "Found clusterings containing %s BGCs across %s cutoffs",
            len(self.computed_labels[thresh]),
            len(self.computed_labels),
        )
