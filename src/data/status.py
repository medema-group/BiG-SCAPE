from src.data.database import Database


def update_bgc_status_all(database: Database, status: int):
    """Updates bgc status"""
    database.update("run_bgc_status",
                    {"status": status})

def update_bgc_status(database: Database, bgc_id: int, status: int):
    """Updates bgc status"""
    database.update("bgc_status",
                    {"status": status},
                    f"where bgc_id = {bgc_id}")

