"""Module that contains helper functions for loading and saving data from and to the cds table in
the database
"""

from src.data.database import Database


def get_cds_rows(database: Database, cds_ids: list = None):
    """Returns a list of rows from the cds table"""
    if cds_ids == None:
        return database.select("cds", "")
    else:
        return database.select("cds", "where id in (" + ",".join(map(str, cds_ids)) + ")")

def gen_accession(base_name, cds_row):
    """generates an accession id for a cds. e.g.
    >AL645882.2.cluster001:gid::pid::loc:12131939:strand:-"""
    nt_start = cds_row["nt_start"]
    nt_end = cds_row["nt_end"]
    strand = cds_row["strand"]
    return f">{base_name}:gid::pid::loc:{nt_start}:{nt_end}:strand:{strand}"