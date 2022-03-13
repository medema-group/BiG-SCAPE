"""Module which contains helper functions for storing/retrieving data related to protein domain
alignments (high scoring protein domains) generated using hmmscan
"""

from src.data.database import Database


def get_predicted_bgc_list(database: Database):
    """Returns a list of ids of genomes which already have predicted bgcs"""
    predicted_bgcs = [
        row["bgc_id"] for row in database.select(
            "bgc_status",
            "where status > 1",
            props=["bgc_id"])]
    return predicted_bgcs

def aa_seq_from_accession(database: Database, accession: str):
    """Returns the amino acid sequences associated with high scoring protein hits by giving
    an accession
    """
    aa_seq = [
        row["aa_seq"] for row in database.select(
            "hsp \
            JOIN hmm ON ",
            "where status > 2",
            props=["cds.aa_seq"])]
    return aa_seq

def insert_hsp(database: Database, cds_id: int, hmm_id: int, bitscore: float):
    """Inserts (or ignores if already there) a new high scoring protein into the database"""
    entry = {
        "cds_id": cds_id,
        "hmm_id": hmm_id,
        "bitscore": bitscore
    }
    database.insert("hsp", entry, True)

def insert_hsp_alignment(database: Database, hsp_id, env_start, env_end, model_start, model_end, model_gaps, cds_start, cds_end, cds_gaps):
    """Inserts (or ignores if already there) a new high scoring protein into the database"""
    entry = {
        "hsp_id": hsp_id,
        "env_start": env_start,
        "env_end": env_end,
        "model_start": model_start,
        "model_end": model_end,
        "model_gaps": model_gaps,
        "cds_start": cds_start,
        "cds_end": cds_end,
        "cds_gaps": cds_gaps
    }
    database.insert("hsp_alignment", entry, True)

def get_hsp_id_list(database):
    """Returns a list of all hsp ids"""
    return [row["id"] for row in database.select("hsp", "", props=["id"])]

def get_hsp_id(database, cds_id, hmm_id):
    """Returns an hsp id based on a given cds id and hmm id"""
    res = database.select("hsp", f"where cds_id = {cds_id} and hmm_id = {hmm_id}", props=["id"])
    if len(res) == 0:
        return None
    else:
        return res[0]["id"]



def get_multiple_align_hsps(database: Database):
    """Returns the rows in the hsp table, joined with the hmm table for accession info
    This also leaves out any rows which have an hmm_id that only occurs once
    """
    # subquery to get a count of hmms (protein domains)
    # done by grouping on hmm_id. if there is only one occurence, it should be removed
    sub_table_query = "hsp \
                       JOIN (SELECT hmm_id AS hmm_count_id, count(*) AS hmmcount \
                       FROM hsp \
                       GROUP BY hmm_count_id)\
                       ON hmm_count_id = hmm_id\
                       JOIN hmm\
                       ON hmm.id = hmm_id"
    rows = database.select(sub_table_query, "WHERE hmmcount > 1", props=["hsp.id", "cds_id", "hmm_id", "hmm.accession"])
    return rows

def generate_pfd_files(run, database: Database):
    """Generates PFS files based on the data in the database"""
    return
