"""Contains tests to test sqlite database reading and writing functions"""

from unittest import TestCase

from sqlalchemy import text

from src.data.sqlite import DB


class TestSQLite(TestCase):
    """Contains tests for sqlite connectivity"""

    def test_create_new_db(self):
        """Tests whether a new db can be created using sqlalchemy"""
        DB.create_new_db()

        get_schema_query = text("SELECT * FROM sqlite_master;")
        result = next(DB.connection.execute(get_schema_query))
        self.assertIsNotNone(result)

        DB.connection.close()

    def test_already_created(self):
        """Tests whether the create database method returns false if a database was
        already created
        """
        DB.create_new_db()

        expected_result = False
        actual_result = DB.create_new_db()

        self.assertEqual(expected_result, actual_result)

    def test_execute_raw_query(self):
        """Tests whether a raw query can be run on a database"""
        DB.create_new_db()

        insert_row_query = (
            "INSERT INTO gbk "
            "(name, as_version, nt_seq, path) "
            "VALUES ('test', 'test', 'test', 'test')"
        )

        result = DB.execute_raw_query(insert_row_query)

        self.assertIsNotNone(result)

        DB.connection.close()
