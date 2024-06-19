# This file contains paths to any included files in the big_scape package

import pkg_resources

DB_SCHEMA_FILE = pkg_resources.resource_filename(
    "bigscape", "big_scape/data/schema.sql"
)
TEMPLATES_OUTPUT_DIR = pkg_resources.resource_filename("bigscape", "big_scape/output")
DEFAULT_CONFIG_FILE = pkg_resources.resource_filename("bigscape", "config.yml")
