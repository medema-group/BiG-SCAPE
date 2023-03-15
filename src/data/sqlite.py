from pathlib import Path
from sqlalchemy import (
    Engine,
    Connection,
    MetaData,
    Compiled,
    CursorResult,
    create_engine,
    text,
)

from src.parameters.constants import DB_SCHEMA_PATH
from src.errors.data import DBClosedError, DBAlreadyOpenError


class DB:
    """Class to manage database interaction"""

    engine: Engine = None
    connection: Connection = None
    metadata: MetaData = None

    @staticmethod
    def opened():
        """Returns true if database is already openened"""
        return DB.connection is not None and not DB.connection.closed

    @staticmethod
    def reflect():
        """Populates the metadata object with information about the tables

        This is necessary since we describe our database in schema.sql and load it, instead
        of creating it through SQLAlchemy's table syste
        """
        if not DB.opened():
            raise DBClosedError()

        DB.metadata = MetaData()
        DB.metadata.reflect(bind=DB.engine)

    @staticmethod
    def create_tables():
        """Populates the database with tables"""
        if not DB.opened:
            raise DBClosedError()

        create_queries = read_schema(DB_SCHEMA_PATH)

        for create_query in create_queries:
            DB.connection.execute(text(create_query))

    @staticmethod
    def create_in_mem():
        """Create a new database in-memory"""

        if DB.opened():
            raise DBAlreadyOpenError()

        DB.engine = create_engine("sqlite:///:memory:")
        DB.connection = DB.engine.connect()

        DB.create_tables()

        DB.reflect()

    @staticmethod
    def close_db():
        """Closes the database connection. This does not save the database to disk"""
        DB.connection.close()

    @staticmethod
    def execute_raw_query(query: str):
        """Executes a raw simple query. Should only be used for very short queries"""
        return DB.connection.execute(text(query))

    @staticmethod
    def execute(query: Compiled) -> CursorResult:
        """Wrapper for SQLAlchemy.connection.execute expecting a Compiled query"""
        return DB.connection.execute(query)


def read_schema(path: Path) -> list[str]:
    with open(path, encoding="utf-8") as schema_file:
        return text_to_queries(schema_file.readlines())


def text_to_queries(schema_lines: list[str]):
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
