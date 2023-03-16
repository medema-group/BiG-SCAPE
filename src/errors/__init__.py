"""Contains custom error classes"""
from src.errors.data import DBAlreadyOpenError, DBClosedError, TableNotFoundError
from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError

__all__ = [
    "DBAlreadyOpenError",
    "DBClosedError",
    "TableNotFoundError",
    "InvalidGBKError",
    "InvalidGBKRegionChildError",
]
