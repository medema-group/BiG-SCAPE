"""Module containing a class for a BiG-SCAPE output object, which has all the
output parameters/arguments"""

# from python
import logging
import time
from typing import Optional
from pathlib import Path

# from other modules
from src.errors import InvalidInputArgError


class OutputParameters:
    """
    Class to store all run output parameters

    Attributes:
        db_path: Path
        log_path: Path
        results_path: Path
    """

    def __init__(self) -> None:
        self.db_path: Optional[Path] = None
        self.log_path: Optional[Path] = None
        self.results_path: Optional[Path] = None

    def parse(self, db_path: Path, log_path: Path, results_path: Path):
        """Load output arguments from commandline ArgParser object

        Args:
            db_path (Path): Path to sqlite database
            log_path (Path): Path to log file
            results_path (Path): Path to output results dir
        """
        self.parse_results_path(results_path)

        # further paths are set to a default path if none are specified
        self.parse_log_path(log_path, results_path)
        self.parse_db_path(db_path, results_path)

    def parse_results_path(self, results_path):
        """Parse the result path parameter and check for errors"""

        if not results_path.parent.exists():
            logging.error("Path to output log file is not valid")
            raise InvalidInputArgError()
        self.results_path = results_path

    def parse_db_path(self, db_path: Path, default_path: Path):
        """Parse the DB path parameter and check for errors

        if no db_path is specified, will set the db_path to a file called data.db
        in the default folder

        if db_path is specified, this function expects that path to point to an
        existing file or a nonexistent location with an existing parent folder. This
        method does not create folders. If db_path points to a folder, this method will
        raise an error
        """

        if db_path is None:
            self.db_path = default_path / Path("data.db")
            return

        # check if parent to a specified path exists
        if db_path and not db_path.parent.exists():
            logging.error("Path to output sqlite db dir is not valid")
            raise InvalidInputArgError()

        if db_path and db_path.is_dir():
            logging.error("Path to output sqlite db dir points to an existing folder")

        self.db_path = db_path

    def parse_log_path(self, log_path: Path, default_path: Path):
        """Parse the log path parameter and check for errors

        If no log_path is specified, will set the log_path to a filename with a
        timestamp created at the time this method is executed, and sets log_path to the
        [default_path]/[timestamp].log

        If log_path is specified, this function expects that path to point to an
        existing file or a nonexistent location with an existing parent folder. This
        method does not create folders. If log_path points to a folder, this method will
        raise an error
        """

        if log_path is None:
            log_timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
            log_filename = log_timestamp + ".log"
            self.log_path = default_path / Path(log_filename)
            return

        # check if parent to a specified path exists
        if log_path and not log_path.parent.exists():
            logging.error("Path to output log file dir is not valid")
            raise InvalidInputArgError()

        if log_path and log_path.is_dir():
            logging.error("Path to output log file points to an existing folder")
            raise InvalidInputArgError()

        self.log_path = log_path
