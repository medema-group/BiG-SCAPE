from sqlalchemy import Engine, Connection, create_engine, text

from src.parameters.constants import DB_SCHEMA_PATH


class DB:
    """Class to manage database interaction"""

    engine: Engine = None
    connection: Connection = None

    @staticmethod
    def create_new_db():
        """Create a new database in-memory"""
        DB.engine = create_engine("sqlite:///:memory:")
        DB.connection = DB.engine.connect()

        with open(DB_SCHEMA_PATH, encoding="utf-8") as schema_file:
            lines = []
            for line in schema_file:
                lines.append(line)
                if not line.rstrip().endswith(";"):
                    continue
                creation_query = text("".join(lines))
                DB.connection.execute(creation_query)
                lines = []

    @staticmethod
    def close_db():
        """Closes the database connection. This does not save the database to disk"""
        DB.connection.close()

    @staticmethod
    def execute_raw_query(query: str):
        """Executes a raw simple query. Should only be used for very short queries"""
        return DB.connection.execute(text(query))
