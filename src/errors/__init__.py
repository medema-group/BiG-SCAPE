"""Contains custom error classes"""
from data import DBAlreadyOpenError, DBClosedError, TableNotFoundError
from genbank import InvalidGBKError, InvalidGBKRegionChildError

__all__ = [
    "DBAlreadyOpenError",
    "DBClosedError",
    "TableNotFoundError",
    "InvalidGBKError",
    "InvalidGBKRegionChildError",
]
