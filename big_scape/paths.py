# This file contains paths to any included files in the big_scape package

import importlib_resources

DB_SCHEMA_FILE = importlib_resources.files("big_scape") / "data/schema.sql"
TEMPLATES_OUTPUT_DIR = importlib_resources.files("big_scape") / "output"
# config is in base directory, not in module
# get path to package directory
DEFAULT_CONFIG_FILE = importlib_resources.files("big_scape") / "config.yml"
