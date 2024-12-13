"""Module containing file loading automation"""
from .load_files import (
    load_dataset_folder,
    get_mibig,
    load_gbks,
    get_all_bgc_records,
    get_all_bgc_records_query,
)

__all__ = [
    "load_dataset_folder",
    "get_mibig",
    "load_gbks",
    "get_all_bgc_records",
    "get_all_bgc_records_query",
]
