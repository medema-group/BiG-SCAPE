"""Contains tests to test sqlite database reading and writing functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from sqlalchemy import text

# from other modules
from big_scape.data import DB
from big_scape.data.sqlite import text_to_queries
from big_scape.errors import DBAlreadyOpenError, DBClosedError


class TestSQLite(TestCase):
    """Contains tests for sqlite connectivity"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_create_in_mem_db(self):
        """Tests whether a new db can be created using sqlalchemy"""
        DB.create_in_mem()

        get_schema_query = text("SELECT * FROM sqlite_master;")
        result = next(DB.connection.execute(get_schema_query))
        self.assertIsNotNone(result)

        DB.close_db()

    def test_reflect_db_not_open(self):
        """Tests whether reflect() raises a DBClosedError when called without
        an active database connection
        """
        self.assertRaises(DBClosedError, DB.reflect)

    def test_create_tables(self):
        """Tests whether create_tables correctly raises a DBClosedError when called
        on a database that is closed or does not exist
        """
        self.assertRaises(DBClosedError, DB.create_tables)

    def test_already_created(self):
        """Tests whether the create database method returns false if a database was
        already created
        """
        DB.open_memory_connection()

        self.assertRaises(DBAlreadyOpenError, DB.create_in_mem)

        DB.close_db()

    def test_execute_raw_query(self):
        """Tests whether a raw query can be run on a database"""
        DB.create_in_mem()

        insert_row_query = (
            "INSERT INTO gbk " "(path, nt_seq) " "VALUES ('test', 'test')"
        )

        result = DB.execute_raw_query(insert_row_query)

        self.assertIsNotNone(result)

        DB.close_db()

    def test_get_table_rows(self):
        """Tests whether the get_table_row_batch function correctly returns a batch
        of rows from the database
        """
        DB.create_in_mem()

        # add 45 rows into database

        for i in range(45):
            insert_row_query = (
                "INSERT INTO gbk " "(path, nt_seq) " f"VALUES ('{i}.gbk', 'test')"
            )
            DB.execute_raw_query(insert_row_query)

        DB.commit()

        expected_row_count = 10

        # get 10 rows from big_scape.database using DB.get_table_rows
        rows = list(next(DB.get_table_row_batch("gbk", expected_row_count)))

        actual_row_count = len(rows)

        self.assertEqual(expected_row_count, actual_row_count)

    def test_save_to_disk(self):
        """Tests whether the sqlite database can be correctly saved to disk"""
        DB.create_in_mem()

        db_save_location = Path("test/test_data/tmp/test.db")

        DB.save_to_disk(db_save_location)

        self.assertTrue(db_save_location.exists())

    def test_save_to_disk_mkdir(self):
        """Tests whether the sqlite database can be correctly saved to disk"""
        DB.create_in_mem()

        db_save_location = Path("test/test_data/tmp/test.db")

        DB.save_to_disk(db_save_location)

        self.assertTrue(db_save_location.parent.exists())

    def test_save_load_from_disk(self):
        """Tests whether the sqlite database can be correctly loaded from disk"""
        db_load_location = Path("test/test_data/database/valid_empty.db")

        DB.load_from_disk(db_load_location)

        expected_tables = ["gbk", "bgc_record"]

        contains_tables = [table in DB.metadata.tables for table in expected_tables]

        self.assertTrue(all(contains_tables))

    def test_text_to_queries(self):
        """Tests the text_to_queries function, which reads a chunk of text and
        returns a list of SQL statements without newlines"""

        # contains 2 queries
        schema_file_path = Path("test/test_data/database/valid_schema.sql")

        with open(schema_file_path, encoding="utf-8") as schema_file:
            query_list = text_to_queries(schema_file.readlines())

        expected_query_count = 2
        actual_query_count = len(query_list)

        self.assertEqual(expected_query_count, actual_query_count)

    def test_reflect(self):
        """Tests whether the metadata correctly reflects the structure of the db
        after calling DB.reflect()
        """
        DB.open_memory_connection()

        # this schema will contain only the gbk and bgc_region tables, and we
        # expect to see them in a moment if the database was correctly reflected
        schema_file_path = Path("test/test_data/database/valid_schema.sql")

        with open(schema_file_path, encoding="utf-8") as schema_file:
            query_list = text_to_queries(schema_file.readlines())

        for query in query_list:
            DB.execute_raw_query(query)

        DB.commit()

        DB.reflect()

        expected_tables = ["gbk", "bgc_record"]

        contains_tables = [table in DB.metadata.tables for table in expected_tables]

        self.assertTrue(all(contains_tables))

        DB.close_db()
