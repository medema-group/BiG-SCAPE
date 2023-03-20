"""Contains custom error classes"""
from src.errors.data import DBAlreadyOpenError, DBClosedError, TableNotFoundError
from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.errors.multithreading import NoFunctionError, WorkerNotStartedError

__all__ = [
    "DBAlreadyOpenError",
    "DBClosedError",
    "TableNotFoundError",
    "InvalidGBKError",
    "InvalidGBKRegionChildError",
    "NoFunctionError",
    "WorkerNotStartedError",
]
