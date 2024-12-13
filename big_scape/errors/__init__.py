"""Contains custom error classes"""
from .data import DBAlreadyOpenError, DBClosedError, TableNotFoundError
from .genbank import InvalidGBKError, InvalidGBKRegionChildError
from .input_args import InvalidArgumentError, ArgumentParseError

__all__ = [
    "DBAlreadyOpenError",
    "DBClosedError",
    "TableNotFoundError",
    "InvalidGBKError",
    "InvalidGBKRegionChildError",
    "InvalidArgumentError",
    "ArgumentParseError",
]
