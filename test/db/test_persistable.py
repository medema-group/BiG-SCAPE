"""Contains tests for the persistable metaclass"""

# from python
from pathlib import Path
from unittest import TestCase

# from other modules
from src.data import Persistable, DB
from src.errors import DBClosedError

# from this module
from .mock_persistable import MockPersistableBGC


class TestPersistable(TestCase):
    """Contains tests for the persistable metaclass"""

    def clean_db(self):
        if DB.opened():
            DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_check_db_connected(self):
        """Tests whether instantiating a new class without starting a new database
        connection throws an error
        """
        self.assertRaises(DBClosedError, Persistable)

    def test_save(self):
        """Tests whether a persistable class saves data to the database correctly"""
        # start db
        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        # MockPersistable is a simple class implementing the persistable class
        mock_persistable = MockPersistableBGC(gbk_file_path)
        mock_persistable.save()

        test_query = "SELECT * FROM gbk;"
        cursor_result = DB.execute_raw_query(test_query)

        gbk_row = cursor_result.mappings().one()

        gbk_row_path = gbk_row["path"]

        self.assertEqual(str(gbk_file_path), gbk_row_path)

        DB.close_db()

    def test_load_one(self):
        """Tests whether a persistable class can be loaded from a database correctly"""

        DB.create_in_mem()

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input.gbk")

        gbk_table = DB.metadata.tables["gbk"]

        insert_query = gbk_table.insert().values(path=str(gbk_file_path)).compile()

        DB.connection.execute(insert_query)

        mock_gbk = MockPersistableBGC.load_one()

        actual_path = mock_gbk.path

        self.assertEqual(gbk_file_path, actual_path)
