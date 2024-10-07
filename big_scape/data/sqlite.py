# from python
from __future__ import annotations
import logging
from datetime import datetime
from pathlib import Path
from typing import Generator, Optional, Any
import sqlite3

# from dependencies
from sqlalchemy import (
    Engine,
    Connection,
    MetaData,
    Compiled,
    Select,
    Insert,
    CursorResult,
    create_engine,
    func,
    select,
    text,
    update,
    desc,
)
from sqlalchemy.pool import StaticPool
import tqdm
import click

# from other modules
import big_scape.paths as bs_paths
from big_scape.cli.config import BigscapeConfig
from big_scape.errors import DBClosedError, DBAlreadyOpenError


class DB:
    """Class to manage database interaction"""

    engine: Optional[Engine] = None
    connection: Optional[Connection] = None
    metadata: Optional[MetaData] = None

    @staticmethod
    def opened() -> bool:
        """Returns true if database is already openened"""
        if DB.engine is None or DB.connection is None:
            return False

        if DB.connection.closed:
            return False

        return True

    @staticmethod
    def reflect() -> None:
        """Populates the metadata object with information about the tables

        This is necessary since we describe our database in schema.sql and load it, instead
        of creating it through SQLAlchemy's table syste
        """
        if not DB.opened():
            raise DBClosedError()

        if DB.engine is None:
            raise RuntimeError("DB.engine is None")

        DB.metadata = MetaData()
        DB.metadata.reflect(bind=DB.engine)

    @staticmethod
    def create_tables() -> None:
        """Populates the database with tables"""
        if not DB.opened():
            raise DBClosedError()

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        creation_queries = read_schema(Path(bs_paths.DB_SCHEMA_FILE))

        for creation_query in creation_queries:
            DB.connection.execute(text(creation_query))

        DB.connection.commit()

    @staticmethod
    def open_memory_connection() -> None:
        """Open a connection to an in-memory database"""

        if DB.opened():
            raise DBAlreadyOpenError()

        DB.engine = create_engine(
            "sqlite:///:memory:",
            connect_args={"check_same_thread": False},
            poolclass=StaticPool,
        )
        DB.connection = DB.engine.connect()

    def open_disk_connection(db_path: Path) -> None:
        if DB.opened():
            raise DBAlreadyOpenError()

        DB.engine = create_engine(
            "sqlite:///" + str(db_path),
            connect_args={"check_same_thread": False},
            poolclass=StaticPool,
        )
        DB.connection = DB.engine.connect()

    @staticmethod
    def create_on_disk(db_path: Path) -> None:
        """Open a connection to a database file"""
        DB.open_disk_connection(db_path)

        DB.create_tables()

        DB.reflect()

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

        # TODO: true arg means silent if there is no context. this is done so that unit
        # tests don't complain. remove this and mock the context in unit tests instead
        click_context = click.get_current_context(silent=True)

        if click_context and click_context.obj["no_db_dump"]:
            return

        # skip this if we are using disk-only mode
        if click_context and click_context.obj["disk_only"]:
            return

        if not DB.opened():
            raise DBClosedError()

        if not DB.engine:
            raise RuntimeError("DB.engine is None")

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        logging.info("Saving database to %s", db_path)

        db_path.parent.mkdir(parents=True, exist_ok=True)

        file_engine = create_engine("sqlite:///" + str(db_path))
        file_engine.connect()

        # from

        raw_memory_connection = DB.engine.raw_connection().driver_connection

        if not isinstance(raw_memory_connection, sqlite3.Connection):
            raise TypeError(
                "Expected raw connection to be of type sqlite3.Connection, got "
                + str(type(raw_memory_connection))
            )

        # to
        raw_file_connection = file_engine.raw_connection().driver_connection

        if not isinstance(raw_file_connection, sqlite3.Connection):
            raise TypeError(
                "Expected raw connection to be of type sqlite3.Connection, got "
                + str(type(raw_file_connection))
            )

        page_count = DB.execute_raw_query("PRAGMA page_count;").scalar_one()

        with tqdm.tqdm(total=page_count) as t:

            def progress(status, remaining, total):
                t.update(total - remaining)

            raw_memory_connection.backup(raw_file_connection, progress=progress)

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

        # disk only means we don't have to dump to memory
        click_context = click.get_current_context(silent=True)
        if click_context and click_context.obj["disk_only"]:
            DB.create_on_disk(db_path)
            return

        file_engine = create_engine("sqlite:///" + str(db_path))
        file_engine.connect()

        DB.open_memory_connection()

        if not DB.engine:
            raise RuntimeError("DB.engine is None")

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        # from
        raw_file_connection = file_engine.raw_connection().driver_connection

        if not isinstance(raw_file_connection, sqlite3.Connection):
            raise TypeError(
                "Expected raw connection to be of type sqlite3.Connection, got "
                + str(type(raw_file_connection))
            )

        # to
        raw_memory_connection = DB.engine.raw_connection().driver_connection

        if not isinstance(raw_memory_connection, sqlite3.Connection):
            raise TypeError(
                "Expected raw connection to be of type sqlite3.Connection, got "
                + str(type(raw_memory_connection))
            )

        # backup only writes those tables that have data, it seems
        DB.create_tables()

        DB.reflect()

        page_count = raw_file_connection.execute("PRAGMA page_count;")
        page_count = page_count.fetchone()[0]

        with tqdm.tqdm(total=page_count, unit="page", desc="Loading database") as t:

            def progress(status, remaining, total):
                t.update(total - remaining)

            raw_file_connection.backup(raw_memory_connection, progress=progress)

    @staticmethod
    def close_db() -> None:
        """Closes the database connection. This does not save the database to disk"""
        if not DB.opened():
            return

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        DB.connection.close()

    @staticmethod
    def execute_raw_query(query: str) -> CursorResult:
        """Executes a raw simple query. Should only be used for very short queries"""
        if not DB.opened():
            raise DBClosedError()

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        return DB.connection.execute(text(query))

    @staticmethod
    def execute(query: Compiled | Select[Any] | Insert, commit=True) -> CursorResult:
        """Wrapper for SQLAlchemy.connection.execute expecting a Compiled query

        Arguments:
            commit: whether or not to immediately commit after executing the query

        This function is meant for single queries.
        """
        if not DB.opened():
            raise DBClosedError()

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        cursor_result = DB.connection.execute(query)  # type: ignore

        if commit:
            DB.connection.commit()

        return cursor_result

    @staticmethod
    def commit() -> None:
        """Performs a commit to the database, saving any alterations to rows and tables
        that have been executed prior

        NOTE: may be redundant if we turn off journaling
        """
        if not DB.opened():
            raise DBClosedError()

        if not DB.connection:
            raise RuntimeError("DB.connection is None")

        DB.connection.commit()

    @staticmethod
    def get_table_row_count(table_name: str) -> int:
        """Return the number of rows in a table

        Args:
            table_name (str): name of the table to get the row count from

        Returns:
            int: number of rows in the table
        """
        if not DB.opened():
            raise DBClosedError()

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

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

        if not DB.opened():
            raise DBClosedError()

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        table_metadata = DB.metadata.tables[table_name]
        table_select = select(table_metadata)
        table_select = table_select.execution_options(stream_results=True)

        cursor = DB.execute(table_select, commit=False)
        while True:
            rows = cursor.fetchmany(batch_size)
            if not rows:
                break
            yield tuple(rows)

    @staticmethod
    def init_run(run: dict) -> None:
        """Initialize a new run in the database

        Args:
            run (dict): run parameters
        """
        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        run_table = DB.metadata.tables["run"]

        run_insert = (
            run_table.insert()
            .values(
                label=run["label"],
                start_time=str(run["start_time"])[:-7],
                mode=run["mode"],
                input_dir=str(run["input_dir"].name) if run["input_dir"] else "None",
                output_dir=str(run["output_dir"].name) if run["output_dir"] else "None",
                reference_dir=(
                    str(run["reference_dir"].name) if run["reference_dir"] else "None"
                ),
                query_path=(
                    str(run["query_bgc_path"].name) if run["query_bgc_path"] else "None"
                ),
                mibig_version=run["mibig_version"] if run["mibig_version"] else "None",
                record_type=run["record_type"].name.title(),
                classify=(
                    "Legacy Groups"
                    if run["legacy_classify"]
                    else (
                        run["classify"].name.title()
                        if run["classify"]
                        else "Not Classify"
                    )
                ),
                weights="Legacy Weights" if run["legacy_weights"] else "Mix",
                alignment_mode=run["alignment_mode"].name.title(),
                extend_strategy=run["extend_strategy"].name.title(),
                include_singletons=(
                    "Yes"
                    if "include_singletons" in run and run["include_singletons"]
                    else "No"
                ),
                cutoffs=",".join(map(str, run["gcf_cutoffs"])),
                config_hash=BigscapeConfig.HASH,
            )
            .returning(run_table.c.id)
            .compile()
        )

        cursor_result = DB.execute(run_insert, False).fetchone()
        if cursor_result is None:
            raise RuntimeError("No run id returned from database!")

        run_id = cursor_result[0]
        run["run_id"] = run_id

    @staticmethod
    def set_run_end(run_id: int, start: datetime, end: datetime):
        """Record run end time and duration to the database

        Args:
            run_id (int): current run id
            start (datetime): run start time
            end (datetime): run end time
        """
        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        run_table = DB.metadata.tables["run"]
        update_stmt = (
            update(run_table)
            .where(run_table.c.id == run_id)
            .values(end_time=str(end)[:-7], duration=str(end - start)[:-7])
        ).compile()
        DB.execute(update_stmt, False)
        DB.commit()

    @staticmethod
    def check_config_hash():
        """Check config file content is the same as the previous run"""
        if DB.metadata is None:
            raise RuntimeError("DB metadata is None")

        run_table = DB.metadata.tables["run"]
        latest_config = DB.execute(
            select(run_table.c.config_hash).order_by(desc(run_table.c.id)).limit(1)
        ).scalar_one()

        if BigscapeConfig.HASH != latest_config:
            raise RuntimeError(
                "Config file values have changed from the previous run! "
                "Existing data is not guarenteed to be reusable, please "
                "run with a fresh output directory/database."
            )


def read_schema(path: Path) -> list[str]:
    """Read an .sql schema from a file"""
    with open(path, encoding="utf-8") as schema_file:
        return text_to_queries(schema_file.readlines())


def text_to_queries(schema_lines: list[str]) -> list[str]:
    """Convert list of lines from an .sql file to a list of queries

    Args:
        schema_lines (list[str]): list of lines from an .sql file

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
