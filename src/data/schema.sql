CREATE TABLE IF NOT EXISTS gbk (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    name VARCHAR,
    as_version TEXT,
    nt_seq TEXT,
    path TEXT,
    UNIQUE(id)
);

-- CREATE TABLE IF NOT EXISTS bgc_region (
--     id INTEGER PRIMARY KEY AUTOINCREMENT,
--     gbk_id INTEGER NOT NULL,
--     number INTEGER,
--     on_contig_edge BOOLEAN,
--     nt_start INTEGER,
--     nt_stop INTEGER,
--     UNIQUE(id),
--     FOREIGN KEY(gbk_id) REFERENCES gbk(id)
-- );
