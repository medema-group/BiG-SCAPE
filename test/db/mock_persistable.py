"""Contains a mock persistable class for testing"""

from src.data import DB, Persistable


class MockPersistableBGC(Persistable):
    def __init__(self, path):
        self.path = path
        return

    def save(self):
        DB.insert_one("gbk", path=str(self.path))
