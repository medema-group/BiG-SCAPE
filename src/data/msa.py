"""Module which contains helper functions for storing/retrieving data related to multi sequence
alignments generated using hmmalign
"""
from src.data.database import Database

def get_aligned_hsp_list(database: Database):
    """Returns a list of all hsp ids in the hsp_alignment table"""
    return [row["hsp_id"] for row in database.select("hsp_alignment", "", props=["hsp_id"])]

def calc_model_gaps():
    return


def calc_cds_gaps():
    return

def insert_msa(database: Database, cds_id, hmm_id, algn_string):
    """Inserts (or ignores if already there) a new entry into the msa table"""
    entry = {
        "cds_id": cds_id,
        "hmm_id": hmm_id,
        "algn_string": algn_string
    }
    database.insert("msa", entry, True)

