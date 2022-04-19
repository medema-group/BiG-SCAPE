import logging
import os
import pyhmmer
from src.big_scape.run.base import Run

from src.data import Database

def get_bigslice_subset():
    """Returns a list of pfam accessions that are included in the bigslice subset of
    accessions, and which will be used later on to do a pre-filter step
    """
    bigslice_accessions = set()
    bigslice_names = set()

    biopfam_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "biopfam.tsv"
    )

    corepfam_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    "corepfam.tsv"
    )

    with open(biopfam_path, encoding="utf-8") as bio_pfam_file:
        for line in bio_pfam_file:
            lineparts = line.rstrip().split("\t")

            if lineparts[3] == "included":
                bigslice_accessions.add(lineparts[0])
                bigslice_names.add(lineparts[1])

    with open(corepfam_path, encoding="utf-8") as core_pfam_file:
        for line in core_pfam_file:
            lineparts = line.rstrip().split("\t")

            bigslice_accessions.add(lineparts[0])
            bigslice_names.add(lineparts[1])

    print(len(bigslice_accessions))

    return bigslice_accessions, bigslice_names

def load_hmms(run, database: Database):
    """Loads the hmm information into the database"""
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    logging.info(" Loading HMM information into database")
    # get hmm profiles

    bigslice_accessions, bigslice_descriptions = get_bigslice_subset()

    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        for profile in hmm_file:
            profile: pyhmmer.plan7.Profile

            bigslice_subset = False

            if profile.accession.decode() in bigslice_accessions:
                bigslice_subset = True

            elif profile.name.decode() in bigslice_descriptions:
                bigslice_subset = True

            database.insert("hmm", {
                "accession": profile.accession.decode(),
                "name": profile.name.decode(),
                "model_length": profile.M
            }, True)
    database.commit_inserts()

def from_id(database: Database, hmm_id: int):
    """Returns a hmm row based on a given ID"""
    return database.select("hmm",
                           f"where id = {hmm_id}",
                           props="*")

def from_accession(database: Database, accession: str):
    """Returns a hmm row based on a given accession"""
    rows = database.select("hmm",
                           f"where accession = \"{accession}\"",
                           props="*")
    if len(rows) == 0:
        return None
    return rows[0]

def gen_subpfam_hmms(run: Run, database: Database):
    return
