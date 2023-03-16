"""Contains an interface class for classes which are committed to the SQLite database"""

# from other modules
from src.errors import DBClosedError

# from this module
from src.data.sqlite import DB


class Persistable:
    """Interface class that describes methods and attributes of a class that persists
    in a database
    """

    def __init__(self):
        Persistable.check_db_connected()

    def save(self):
        """Saves the current object to a database"""
        raise NotImplementedError()

    @classmethod
    def load_one(cls):
        """Loads an instance of the current object from a database"""
        raise NotImplementedError()

    @classmethod
    def load_all(cls):
        """Loads all instances found in the database"""
        raise NotImplementedError()

    @staticmethod
    def check_db_connected():
        if not DB.opened():
            raise DBClosedError()
