"""Contains tests to test sqlite database reading and writing functions"""

from unittest import TestCase

from src.data.db import DB


class TestSQLite(TestCase):
    """Contains tests for sqlite connectivity"""

    def test_create_new_db(self):
        """Tests whether a new db can be created using sqlalchemy"""
        DB.create_new_db()

        DB.execute_raw_query(
            "INSERT INTO gbk"
            "(name, as_version, nt_seq, path)"
            "VALUES ('test', 'test', 'test', 'test');"
        )

        DB.execute_raw_query("SELECT * FROM gbk;")
