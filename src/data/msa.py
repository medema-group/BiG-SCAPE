"""Module which contains helper functions for storing/retrieving data related
to multi sequence alignments generated using hmmalign

Author: Arjan Draisma
"""
from src.data.database import Database

def get_aligned_hsp_list(database: Database):
    """Returns a list of all hsp ids in the hsp_alignment table"""
    return [row["id"]
        for row in database.select(
            "hsp \
            join msa on msa.cds_id = hsp.cds_id \
            and msa.hmm_id = hsp.hmm_id",
            "",
            props=["hsp.id"]
    )]

def insert_msa(database: Database, cds_id, hmm_id, env_start, env_end, algn_string):
    """Inserts (or ignores if already there) a new entry into the msa table"""
    entry = {
        "cds_id": cds_id,
        "hmm_id": hmm_id,
        "env_start": env_start,
        "env_end": env_end,
        "algn_string": algn_string
    }
    database.insert("msa", entry, True)

