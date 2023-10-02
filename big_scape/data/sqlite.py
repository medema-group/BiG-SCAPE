# from python
import logging
from pathlib import Path
from typing import Generator, List
from sqlite3 import Connection as Sqlite3Connection

# from dependencies
from sqlalchemy import (
    Engine,
    Connection,
    MetaData,
    Compiled,
    CursorResult,
    create_engine,
    func,
    select,
    text,
)

# from other modules
from big_scape.parameters.constants import DB_SCHEMA_PATH
from big_scape.errors import DBClosedError, DBAlreadyOpenError


class DB:
    """Class to manage database interaction"""

    engine: Engine = None
    connection: Connection = None
    metadata: MetaData = None

    @staticmethod
    def opened() -> bool:
        """Returns true if database is already openened"""
        return DB.connection is not None and not DB.connection.closed

    @staticmethod
    def reflect() -> None:
        """Populates the metadata object with information about the tables

        This is necessary since we describe our database in schema.sql and load it, instead
        of creating it through SQLAlchemy's table syste
        """
        if not DB.opened():
            raise DBClosedError()

        DB.metadata = MetaData()
        DB.metadata.reflect(bind=DB.engine)

    @staticmethod
    def create_tables() -> None:
        """Populates the database with tables"""
        if not DB.opened():
            raise DBClosedError()

        creation_queries = read_schema(Path(DB_SCHEMA_PATH))

        for creation_query in creation_queries:
            DB.connection.execute(text(creation_query))

        DB.connection.commit()

    @staticmethod
    def open_memory_connection() -> None:
        """Open a connection to an in-memory database"""

        if DB.opened():
            raise DBAlreadyOpenError()

        DB.engine = create_engine("sqlite:///:memory:")
        DB.connection = DB.engine.connect()

    @staticmethod
    def create_in_mem() -> None:
        """Create a new database in-memory"""
        DB.open_memory_connection()

        DB.create_tables()

        DB.reflect()

    @staticmethod
    def save_to_disk(db_path: Path) -> None:
        """Saves the in-memory database to a .db file. This overwrites any last database
        file in the same location
        """
        if not DB.opened():
            raise DBClosedError()

        logging.info("Saving database to %s", db_path)

        db_path.parent.mkdir(parents=True, exist_ok=True)

        file_engine = create_engine("sqlite:///" + str(db_path))
        file_engine.connect()

        # from
        raw_memory_connection: Sqlite3Connection = (
            DB.engine.raw_connection().driver_connection
        )
        # to
        raw_file_connection: Sqlite3Connection = (
            file_engine.raw_connection().driver_connection
        )

        raw_memory_connection.backup(raw_file_connection)

    @staticmethod
    def load_from_disk(db_path: Path) -> None:
        """Loads the database from a database file to memory

        Args:
            db_path (Path): Path to the database file to load
        """
        if DB.opened():
            raise DBAlreadyOpenError()

        if not db_path.exists():
            raise FileNotFoundError()

        file_engine = create_engine("sqlite:///" + str(db_path))
        file_engine.connect()

        DB.open_memory_connection()

        # from
        raw_file_connection: Sqlite3Connection = (
            file_engine.raw_connection().driver_connection
        )
        # to
        raw_memory_connection: Sqlite3Connection = (
            DB.engine.raw_connection().driver_connection
        )

        # backup only writes those tables that have data, it seems
        DB.create_tables()

        DB.reflect()

        raw_file_connection.backup(raw_memory_connection)

    @staticmethod
    def close_db() -> None:
        """Closes the database connection. This does not save the database to disk"""
        DB.connection.close()

    @staticmethod
    def execute_raw_query(query: str) -> CursorResult:
        """Executes a raw simple query. Should only be used for very short queries"""
        return DB.connection.execute(text(query))

    @staticmethod
    def execute(query: Compiled, commit=True) -> CursorResult:
        """Wrapper for SQLAlchemy.connection.execute expecting a Compiled query

        Arguments:
            commit: whether or not to immediately commit after executing the query

        This function is meant for single queries.
        """

        cursor_result = DB.connection.execute(query)

        if commit:
            DB.connection.commit()

        return cursor_result

    @staticmethod
    def commit() -> None:
        """Performs a commit to the database, saving any alterations to rows and tables
        that have been executed prior

        NOTE: may be redundant if we turn off journaling
        """
        DB.connection.commit()

    @staticmethod
    def get_table_row_count(table_name: str) -> int:
        """Return the number of rows in a table

        Args:
            table_name (str): name of the table to get the row count from

        Returns:
            int: number of rows in the table
        """
        table_metadata = DB.metadata.tables[table_name]
        return DB.execute(
            select(func.count("*")).select_from(table_metadata)
        ).scalar_one()

    @staticmethod
    def get_table_row_batch(
        table_name: str, batch_size=100
    ) -> Generator[tuple, None, None]:
        """Generator that yields rows from a table

        Args:
            table_name (str): name of the table to get rows from
            batch_size (int): how many rows to yield at a time

        Yields:
            Generator[tuple]: generator of tuples of rows from the table
        """
        table_metadata = DB.metadata.tables[table_name]
        table_select = select(table_metadata)
        table_select = table_select.execution_options(stream_results=True)

        cursor = DB.execute(table_select, commit=False)
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                break
            yield rows


def read_schema(path: Path) -> List[str]:
    """Read an .sql schema from a file"""
    with open(path, encoding="utf-8") as schema_file:
        return text_to_queries(schema_file.readlines())


def text_to_queries(schema_lines: List[str]) -> list[str]:
    """Convert list of lines from an .sql file to a list of queries

    Args:
        schema_lines (List[str]): list of lines from an .sql file

    Returns:
        list[str]: list of queries that can be executed
    """
    create_queries = []
    query_lines = []
    for line in schema_lines:
        query_lines.append(line)
        if not line.rstrip().endswith(";"):
            continue
        creation_query = "".join(query_lines)
        create_queries.append(creation_query)
        query_lines = []
    return create_queries
