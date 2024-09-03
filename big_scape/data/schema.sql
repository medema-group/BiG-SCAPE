-- pragmas


-- regular tables
CREATE TABLE IF NOT EXISTS run (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    label TEXT,
    start_time TEXT,
    end_time TEXT,
    duration TEXT,
    mode TEXT, --cluster/query
    input_dir TEXT,
    output_dir TEXT,
    reference_dir TEXT,
    query_path TEXT,
    mibig_version TEXT,
    record_type TEXT,
    classify TEXT,
    weights TEXT,
    alignment_mode TEXT,
    include_singletons TEXT,
    cutoffs TEXT
);

CREATE TABLE IF NOT EXISTS gbk (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    path TEXT,
    hash TEXT,
    nt_seq TEXT,
    organism TEXT,
    taxonomy TEXT,
    description TEXT,
    UNIQUE(hash)
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
    merged BOOLEAN,
    UNIQUE(id),
    UNIQUE(gbk_id, record_number, record_type),
    FOREIGN KEY(gbk_id) REFERENCES gbk(id)
);

CREATE TABLE IF NOT EXISTS bgc_record_family (
    record_id INTEGER NOT NULL,
    family_id INTEGER NOT NULL,
    UNIQUE(record_id, family_id),
    FOREIGN KEY(record_id) REFERENCES bgc_record(id)
    FOREIGN KEY(family_id) REFERENCES family(id)
);

CREATE TABLE IF NOT EXISTS family (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    center_id INTEGER NOT NULL,
    newick TEXT,
    cutoff REAL NOT NULL,
    bin_label TEXT NOT NULL,
    run_id INTEGER NOT NULL,
    UNIQUE(id),
    UNIQUE(center_id, cutoff, bin_label, run_id),
    FOREIGN KEY(run_id) REFERENCES run(id)
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
    record_a_id INTEGER NOT NULL,
    record_b_id INTEGER NOT NULL,
    distance REAL NOT NULL,
    jaccard REAL NOT NULL,
    adjacency REAL NOT NULL,
    dss REAL NOT NULL,
    edge_param_id INTEGER NOT NULL,
    lcs_a_start INTEGER NOT NULL,
    lcs_a_stop INTEGER NOT NULL,
    lcs_b_start INTEGER NOT NULL,
    lcs_b_stop INTEGER NOT NULL,
    ext_a_start INTEGER NOT NULL,
    ext_a_stop INTEGER NOT NULL,
    ext_b_start INTEGER NOT NULL,
    ext_b_stop INTEGER NOT NULL,
    reverse BOOLEAN NOT NULL,
    lcs_domain_a_start INTEGER NOT NULL,
    lcs_domain_a_stop INTEGER NOT NULL,
    lcs_domain_b_start INTEGER NOT NULL,
    lcs_domain_b_stop INTEGER NOT NULL,
    UNIQUE(record_a_id, record_b_id, edge_param_id)
    FOREIGN KEY(record_a_id) REFERENCES bgc_record(id)
    FOREIGN KEY(record_b_id) REFERENCES bgc_record(id)
    FOREIGN KEY(edge_param_id) REFERENCES edge_params(id)
);

CREATE TABLE IF NOT EXISTS connected_component (
    id INTEGER NOT NULL,
    record_id INTEGER NOT NULL,
    cutoff REAL NOT NULL,
    edge_param_id INTEGER NOT NULL,
    bin_label TEXT NOT NULL,
    run_id INTEGER NOT NULL,
    UNIQUE(record_id, cutoff, bin_label, run_id)
    FOREIGN KEY(record_id) REFERENCES bgc_record(id)
    FOREIGN KEY(run_id) REFERENCES run(id)
);

CREATE TABLE IF NOT EXISTS edge_params (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    weights TEXT NOT NULL,
    alignment_mode TEXT NOT NULL,
    UNIQUE(id),
    UNIQUE(weights, alignment_mode)
);

CREATE INDEX IF NOT EXISTS record_id_index ON bgc_record(id);
CREATE INDEX IF NOT EXISTS distance_record_id_index ON distance(record_a_id);
CREATE INDEX IF NOT EXISTS distance_record_id_index ON distance(record_b_id);
