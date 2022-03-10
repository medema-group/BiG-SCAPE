import logging
import os
import pyhmmer

from src.data.database import Database


def load_hmms(run, database: Database):
    """Loads the hmm information into the database"""
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    logging.info("Loading HMM information into database")
    # get hmm profiles
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
    
        for profile in hmm_file:
            profile: pyhmmer.plan7.Profile
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
