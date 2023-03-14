"""Contains errors raised for data handling using SQLite"""


class DatabaseNotCreatedError(Exception):
    def __init__(self):
        super().__init__(
            (
                "Database connection not created before operation."
                "Create a new database or load it from disk first"
            )
        )


class DatabaseClosedError(Exception):
    def __init__(self):
        super().__init__(
            (
                "Database was created before, but is closed now."
                "Did you close the database in your code before trying to read/write it?"
            )
        )
