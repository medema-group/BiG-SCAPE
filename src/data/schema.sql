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

CREATE TABLE IF NOT EXISTS bgc_region (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    parent_id INTEGER NOT NULL,
    region_number INTEGER,
    gbk_id INTEGER NOT NULL,
    on_contig_edge BOOLEAN,
    nt_start INTEGER,
    nt_stop INTEGER,
    UNIQUE(id),
    UNIQUE(parent_id, region_number)
    FOREIGN KEY(parent_id) REFERENCES bgc_region(id),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS bgc_region_type (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    region_id INTEGER NOT NULL,
    type TEXT NOT NULL,
    UNIQUE(region_id, type),
    FORIEGN KEY (region_id) REFERENCES bgc_region(id)
)

CREATE TABLE IF NOT EXISTS cds (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    region_id INTEGER NOT NULL,
    nt_start INTEGER,
    nt_stop INTEGER,
    strand INTEGER,
    locus_tag TEXT,
    protein_id TEXT,
    product TEXT,
    aa_seq TEXT,
    UNIQUE(id, region_id),
    FOREIGN KEY(region_id) REFERENCES bgc_region(id)
)
