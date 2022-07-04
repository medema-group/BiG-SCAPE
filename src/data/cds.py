"""Module that contains helper functions for loading and saving data from and
to the cds table in a database

Authors: Satria A. Kautsar, Arjan Draisma.

Copied from
https://github.com/medema-group/bigslice
file: bigslice/modules/data/bgc.py

Modified by Arjan Draisma to work with BiG-SCAPE
"""

from src.data.database import Database


def get_cds_rows(database: Database, cds_ids: list = None):
    """Returns a list of rows from the cds table"""
    if cds_ids == None:
        return database.select("cds", "")
    else:
        return database.select("cds", "where id in (" + ",".join(map(str, cds_ids)) + ")")

def get_cds_with_alignment(database: Database, bgc_name = None):
    """Returns a complete list of cds entries and protein domain alignment details"""
    if bgc_name is None:
        clause = ""
    else:
        clause = f"where bgc.name = \"{bgc_name}\""
    return database.select(
        "hsp_alignment \
        join hsp on hsp.id = hsp_alignment.hsp_id \
        join cds on cds.id = hsp.cds_id \
        join bgc on bgc.id = cds.bgc_id \
        join hmm on hmm.id = hsp.hmm_id",
        clause,
        props=[
            "bgc.name",
            "cds.orf_id",
            "cds.nt_start",
            "cds.nt_end",
            "cds.strand",
            "hmm.accession",
            "hsp_alignment.env_start",
            "hsp_alignment.env_end",
            "hsp.bitscore"
        ]
    )

def gen_header(base_name, cds_row):
    """generates an accession id for a cds. e.g.
    >AL645882.2.cluster001:gid::pid::loc:12131939:strand:-
    From a cds row returned from the database
    """
    nt_start = cds_row["nt_start"]
    nt_end = cds_row["nt_end"]
    # TODO: replace 1 for forward and -1 for backwards with + and -
    strand = cds_row["strand"]
    return f">{base_name}:gid::pid::loc:{nt_start}:{nt_end}:strand:{strand}"

def gen_header_cds(base_name, cds_obj):
    """generates an accession id for a cds. e.g.
    >AL645882.2.cluster001:gid::pid::loc:12131939:strand:-
    from a cds object
    """
    nt_start = cds_obj.nt_start
    nt_end = cds_obj.nt_end
    strand = "+" if cds_obj.strand == 1 else "-"
    return f">{base_name}:gid::pid::loc:{nt_start}:{nt_end}:strand:{strand}"


def get_aa_from_header(database: Database, header):
    """Gets the amino acid sequence from a header"""
    parts = header.split(":")
    name = parts[0]
    start = parts[6]
    end = parts[7]
    strand = parts[9]
    rows = database.select(
        "cds \
        join bgc on bgc.id = cds.bgc_id",
        f"where name = \"{name}\" \
        and nt_start = {start} \
        and nt_end = {end} \
        and strand = {strand}",
        props=["aa_seq"]
    )
    if len(rows) == 0:
        return None
    else:
        return rows[0]["aa_seq"]
