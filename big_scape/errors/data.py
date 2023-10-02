"""Contains errors raised for data handling using SQLite"""


class DBClosedError(Exception):
    """Thrown when the database is closed. This means the table was closed after creation
    or was never created in the first place
    """

    def __init__(self):
        super().__init__(
            (
                "Database was created before, but is closed now."
                "Did you close the database in your code before trying to read/write it?"
            )
        )


class DBAlreadyOpenError(Exception):
    """Error thrown when trying to open a database when a database is already open"""

    def __init__(self):
        super().__init__("Database already opened!")


class TableNotFoundError(Exception):
    """Error thrown in any case that a table name is given to a function that does not exist"""

    def __init__(self, table_name):
        super().__init__(f"Table {table_name} does not exist in database!")
