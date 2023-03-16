"""Contains tests to test sqlite database reading and writing functions"""

# from python
from pathlib import Path
from unittest import TestCase

# from dependencies
from sqlalchemy import text

# from other modules
from src.data import DB
from src.errors import DBAlreadyOpenError


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
