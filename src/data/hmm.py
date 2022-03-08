


def (bgc_ids: List[int], database: Database):
    """query database, get all aa sequences
    of the CDS into a multifasta string
    e.g. for the purpose of doing hmmscan"""

    rows = database.select(
        "cds",
        "WHERE bgc_id IN (" + ",".join(map(str, bgc_ids)) + ")",
        props=["nt_start", "nt_end", "strand", "locus_tag", "protein_id", "product", "aa_seq"]
    )

    return rows