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
    parent_id INTEGER NOT NULL,
    record_number INTEGER,
    gbk_id INTEGER NOT NULL,
    on_contig_edge BOOLEAN,
    nt_start INTEGER,
    nt_stop INTEGER,
    UNIQUE(id),
    UNIQUE(parent_id, record_number)
    FOREIGN KEY(parent_id) REFERENCES bgc_record(id),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id) );

-- another comment
