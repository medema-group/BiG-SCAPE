"""Contains errors raised for data handling using SQLite"""

from errors import Exception


class DatabaseNotCreatedError(Exception):
    def __init__(self):
        super().__init__(
            (
                "Database connection not created before operation."
                "Create a new database or load it from disk first"
            )
        )
