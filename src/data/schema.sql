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
    record_number INTEGER,
    gbk_id INTEGER,
    contig_edge BOOLEAN,
    type TEXT,
    nt_start INTEGER,
    nt_stop INTEGER,
    UNIQUE(id),
    UNIQUE(gbk_id, type, nt_start, nt_stop),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS cds (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gbk_id INTEGER NOT NULL,
    nt_start INTEGER NOT NULL,
    nt_stop INTEGER NOT NULL,
    orf_num INTEGER NOT NULL,
    strand INTEGER NOT NULL,
    gene_kind TEXT,
    aa_seq TEXT NOT NULL,
    UNIQUE(id),
    UNIQUE(gbk_id, nt_start, nt_stop, orf_num, strand)
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS hsp (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    cds_id INTEGER NOT NULL,
    hmm_id INTEGER NOT NULL,
    env_start INTEGER NOT NULL,
    env_stop INTEGER NOT NULL,
    bitscore REAL NOT NULL
);

CREATE TABLE IF NOT EXISTS hsp_alignment (
    hsp_id INTEGER PRIMARY KEY NOT NULL,
    alignment TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS hmm (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    accession TEXT NOT NULL,
    name TEXT NOT NULL,
    db_id INTEGER NOT NULL,
    model_length INTEGER NOT NULL,
    FOREIGN KEY(db_id) REFERENCES hmm_db(id)
);

CREATE TABLE IF NOT EXISTS hmm_db (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    path TEXT NOT NULL,
    md5 TEXT NOT NULL
);
