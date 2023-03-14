"""Contains tests to test sqlite database reading and writing functions"""

from pathlib import Path
from unittest import TestCase

from sqlalchemy import text

from src.data.sqlite import DB
from src.errors.data import DBAlreadyOpenError, TableNotFoundError


class TestSQLite(TestCase):
    """Contains tests for sqlite connectivity"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_read_schema(self):
        """Tests whether a schema detailing table creation is correctly read"""
        pass

    def test_create_in_mem_db(self):
        """Tests whether a new db can be created using sqlalchemy"""
        DB.create_in_mem()

        get_schema_query = text("SELECT * FROM sqlite_master;")
        result = next(DB.connection.execute(get_schema_query))
        self.assertIsNotNone(result)

        DB.close_db()

    def test_already_created(self):
        """Tests whether the create database method returns false if a database was
        already created
        """
        DB.create_in_mem()

        self.assertRaises(DBAlreadyOpenError, DB.create_in_mem)

        DB.close_db()

    def test_execute_raw_query(self):
        """Tests whether a raw query can be run on a database"""
        DB.create_in_mem()

        insert_row_query = (
            "INSERT INTO gbk "
            "(name, as_version, nt_seq, path) "
            "VALUES ('test', 'test', 'test', 'test')"
        )

        result = DB.execute_raw_query(insert_row_query)

        self.assertIsNotNone(result)

        DB.close_db()

    def test_insert_one(self):
        """Tests the insert_one function, specifically meant to immediately insert
        one row into a given table
        """
        DB.create_in_mem()

        table_name = "gbk"
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        DB.insert_one(table_name, path=str(gbk_file_path))

        cursor_result = DB.connection.execute(text("SELECT * FROM gbk;"))

        expected_row_count = 1
        actual_row_count = len(cursor_result.fetchall())

        self.assertEqual(expected_row_count, actual_row_count)

        DB.close_db()

    def test_insert_one_table_does_not_exist(self):
        """Tests whether insert_one correctly throws an error when a table is given that
        does not exist
        """
        DB.create_in_mem()

        target_table = "fakename"
        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        self.assertRaises(
            TableNotFoundError, DB.insert_one, target_table, path=gbk_file_path
        )

        DB.close_db()
