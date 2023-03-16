"""Contains a mock persistable class for testing"""

from pathlib import Path
from src.data import DB, Persistable


class MockPersistableBGC(Persistable):
    def __init__(self, path):
        self.path = Path(path)
        return

    def save(self):
        gbk_table = DB.metadata.tables["gbk"]

        insert_query = gbk_table.insert().values(path=str(self.path)).compile()

        DB.connection.execute(insert_query)

    @classmethod
    def load_one(cls):
        gbk_table = DB.metadata.tables["gbk"]

        select_query = gbk_table.select().limit(1).compile()

        cursor_result = DB.execute(select_query)

        gbk_map = cursor_result.mappings().one()

        gbk_path = gbk_map["path"]

        new_gbk = cls(gbk_path)

        return new_gbk
