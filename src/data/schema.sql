-- pragmas


-- regular tables
CREATE TABLE IF NOT EXISTS gbk (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name VARCHAR,
    as_version TEXT,
    nt_seq TEXT,
    path TEXT,
    UNIQUE(id)
);

CREATE TABLE IF NOT EXISTS bgc_record (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    parent_id INTEGER,
    record_number INTEGER,
    gbk_id INTEGER,
    contig_edge BOOLEAN,
    type TEXT,
    nt_start INTEGER,
    nt_stop INTEGER,
    UNIQUE(id),
    UNIQUE(parent_id, record_number)
    FOREIGN KEY(parent_id) REFERENCES bgc_record(id),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS cds (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    record_id INTEGER NOT NULL,
    nt_start INTEGER NOT NULL,
    nt_stop INTEGER NOT NULL,
    strand INTEGER NOT NULL,
    locus_tag TEXT NOT NULL,
    protein_id TEXT,
    product TEXT,
    aa_seq TEXT NOT NULL,
    UNIQUE(id, record_id),
    FOREIGN KEY(record_id) REFERENCES bgc_record(id)
)

CREATE TABLE IF NOT EXISTS hsp (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    cds_id INTEGER NOT NULL,
    hmm_id INTEGER NOT NULL,
    bitscore REAL NOT NULL
)

CREATE TABLE IF NOT EXISTS hsp_alignment (
    hsp_id INTEGER PRIMARY KEY NOT NULL,
    model_start INTEGER NOT NULL,
    model_stop INTEGER NOT NULL,
    model_gaps TEXT NOT NULL,
    cds_start INTEGER NOT NULL,
    cds_stop INTEGER NOT NULL,
    cds_gaps TEXT NOT NULL
)

CREATE TABLE IF NOT EXISTS hmm (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    accession TEXT NOT NULL,
    name TEXT NOT NULL,
    db_id INTEGER NOT NULL,
    model_length INTEGER NOT NULL,
    FOREIGN KEY(db_id) REFERENCES hmm_db(id)
)

CREATE TABLE IF NOT EXISTS hmm_db (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    path TEXT NOT NULL,
    md5 TEXT NOT NULL
)
