"""Contains custom error classes"""
from src.errors.data import DBAlreadyOpenError, DBClosedError, TableNotFoundError
from src.errors.genbank import InvalidGBKError, InvalidGBKRegionChildError
from src.errors.multithreading import (
    WorkerPoolSetupError,
    WorkerSetupError,
    WorkerExecutionError,
)
from src.errors.input_args import InvalidInputArgError

__all__ = [
    "DBAlreadyOpenError",
    "DBClosedError",
    "TableNotFoundError",
    "InvalidGBKError",
    "InvalidGBKRegionChildError",
    "WorkerPoolSetupError",
    "WorkerSetupError",
    "WorkerExecutionError",
    "InvalidInputArgError",
]
