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
