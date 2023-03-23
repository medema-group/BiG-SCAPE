"""Module containing a class for a BiG-SCAPE output object, which has all the
output parameters/arguments"""

# from python
import logging
from typing import Optional
from pathlib import Path

# from dependencies
# from other modules
from src.errors import InvalidInputArgError

# from this module


class Output:
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

        if not results_path.parent.exists():
            logging.error("Path to output log file is not valid")
            raise InvalidInputArgError()
        self.results_path = results_path

        if not db_path:
            self.db_path = results_path

        if db_path and not db_path.parent.exists():
            logging.error("Path to output sqlite db dir is not valid")
            raise InvalidInputArgError()
        self.db_path = db_path

        if not log_path:
            self.log_path = results_path

        if log_path and not log_path.parent.exists():
            logging.error("Path to output log file dir is not valid")
            raise InvalidInputArgError()
        self.log_path = log_path
