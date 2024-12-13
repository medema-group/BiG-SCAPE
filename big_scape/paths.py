# This file contains paths to any included files in the big_scape package

import pkg_resources

DB_SCHEMA_FILE = pkg_resources.resource_filename("big_scape", "data/schema.sql")
TEMPLATES_OUTPUT_DIR = pkg_resources.resource_filename("big_scape", "output")
# config is in base directory, not in module
# get path to package directory
DEFAULT_CONFIG_FILE = pkg_resources.resource_filename("big_scape", "config.yml")
