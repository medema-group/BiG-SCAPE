-- pragmas


-- regular tables
CREATE TABLE IF NOT EXISTS gbk (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    path TEXT,
    source_type TEXT,
    nt_seq TEXT,
    UNIQUE(path)
);

CREATE TABLE IF NOT EXISTS bgc_record (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    gbk_id INTEGER,
    parent_id INTEGER,
    record_number INTEGER,
    contig_edge BOOLEAN,
    record_type TEXT,
    nt_start INTEGER,
    nt_stop INTEGER,
    product TEXT,
    category TEXT,
    UNIQUE(id),
    UNIQUE(gbk_id, record_number, record_type),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS bgc_record_family (
    record_id INTEGER NOT NULL,
    family INTEGER NOT NULL,
    cutoff REAL NOT NULL,
    UNIQUE(record_id, cutoff, family),
    FOREIGN KEY(record_id) REFERENCES bgc_record(id)
);

CREATE TABLE IF NOT EXISTS scanned_cds (
    cds_id INTEGER PRIMARY KEY NOT NULL,
    unique(cds_id),
    FOREIGN KEY(cds_id) REFERENCES cds(id)
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
    accession TEXT NOT NULL,
    env_start INTEGER NOT NULL,
    env_stop INTEGER NOT NULL,
    bit_score REAL NOT NULL
);

CREATE TABLE IF NOT EXISTS hsp_alignment (
    hsp_id INTEGER PRIMARY KEY NOT NULL,
    alignment TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS distance (
    region_a_id INTEGER NOT NULL,
    region_b_id INTEGER NOT NULL,
    distance REAL NOT NULL,
    jaccard REAL NOT NULL,
    adjacency REAL NOT NULL,
    dss REAL NOT NULL,
    weights TEXT NOT NULL,
    UNIQUE(region_a_id, region_b_id)
    FOREIGN KEY(region_a_id) REFERENCES bgc_record(id)
    FOREIGN KEY(region_b_id) REFERENCES bgc_record(id)
);
