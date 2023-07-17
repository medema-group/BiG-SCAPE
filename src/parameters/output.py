"""Module containing a class for a BiG-SCAPE output object, which has all the
output parameters/arguments"""

# from python
import logging
import time
from pathlib import Path

# from other modules
from src.errors import InvalidArgumentError


class OutputParameters:
    """Class to store all output parameters

    Attributes:
        output_dir: Path
        db_path: Path
        log_path: Path
        profile_path: Path

    """

    def __init__(self):
        self.output_dir: Path = Path("")
        self.db_path: Path = Path("")
        self.log_path: Path = Path("")
        self.profile_path: Path = Path("")

    def validate(self):
        """Load output arguments from commandline ArgParser object

        Args:
            db_path (Path): Path to sqlite database
            log_path (Path): Path to log file
            results_path (Path): Path to output results dir
        """
        # generate timestamp
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

        validate_output_dir(self.output_dir)

        # further paths are set to a default path if none are specified
        self.db_path = validate_db_path(self.db_path, self.output_dir)
        self.log_path = validate_log_path(self.log_path, self.output_dir, timestamp)
        self.profile_path = validate_profile_path(
            self.profile_path, self.output_dir, timestamp
        )


def validate_output_dir(output_dir: Path) -> None:
    """Parse the result path parameter and check for errors"""

    # parent must exist
    if not output_dir.parent.exists():
        logging.error("Path to output log file is not valid")
        raise InvalidArgumentError("--output_dir", output_dir)

    # create output directory
    if output_dir.exists():
        logging.info("Reusing output directory: %s", output_dir)
        return

    logging.info("Creating output directory: %s", output_dir)
    output_dir.mkdir()


def validate_db_path(db_path: Path, default_path: Path) -> Path:
    """Parse the DB path parameter and check for errors

    if no db_path is specified, will set the db_path to a file called data.db
    in the default folder in the format:
    [default_path]/data.db

    if db_path is specified, this function expects that path to point to an
    existing file or a nonexistent location with an existing parent folder. This
    method does not create folders. If db_path points to a folder, this method will
    raise an error
    """

    # set to default if this was initially None
    if db_path is None:
        # but we do need a default path. This exception should not happen since we are
        # validating the default path before this function, and it should never be None
        # anyway
        if default_path is None:
            logging.error("Default path is not set!")
            raise InvalidArgumentError("--db_path", db_path)

        return default_path / Path("data_sqlite.db")

    # check if parent to a specified path exists
    if not db_path.parent.exists():
        logging.error("Path to output sqlite db dir is not valid")
        raise InvalidArgumentError("--db_path", db_path)

    # db_path should point to an existing or empty file. if a folder exists at the
    # path instead we throw an error
    if db_path.is_dir():
        logging.error("Path to output sqlite db dir points to an existing folder")
        raise InvalidArgumentError("--db_path", db_path)

    return db_path


def validate_log_path(log_path: Path, default_path: Path, timestamp: str) -> Path:
    """Validate the log_path parameter

    If no log_path is specified, will return a path with a filename consisting of a
    timestamp created at the time this method is executed, in the format:
    [default_path]/[timestamp].log

    If log_path is specified, this function expects that path to point to an
    existing file or a nonexistent location with an existing parent folder. This
    method does not create folders. If log_path points to a folder, this method will
    raise an error
    """

    if log_path is None:
        # but we do need a default path. This exception should not happen since we are
        # validating the default path before this function, and it should never be None
        # anyway
        if default_path is None:
            logging.error("Default path is not set!")
            raise InvalidArgumentError("--log_path", log_path)

        # return new path for log
        log_filename = timestamp + ".log"
        return default_path / Path(log_filename)

    # check if parent to a specified path exists
    if log_path and not log_path.parent.exists():
        logging.error("Path to output log file dir is not valid")
        raise InvalidArgumentError("--log_path", log_path)

    if log_path and log_path.is_dir():
        logging.error("Path to output log file points to an existing folder")
        raise InvalidArgumentError("--log_path", log_path)

    # otherwise we're good!
    return log_path


def validate_profile_path(
    profile_path: Path, default_path: Path, timestamp: str
) -> Path:
    """Validate the profile_path parameter

    If no profile_path is specified, will return a path with a filename consisting of a
    timestamp created at the time this method is executed, in the format:
    [default_path]/[timestamp].profile

    If profile_path is specified, this function expects that path to point to an
    existing file or a nonexistent location with an existing parent folder. This
    method does not create folders. If profile_path points to a folder, this method will
    raise an error
    """

    if profile_path is None:
        # but we do need a default path. This exception should not happen since we are
        # validating the default path before this function, and it should never be None
        # anyway
        if default_path is None:
            logging.error("Default path is not set!")
            raise InvalidArgumentError("--log_path", profile_path)

        # return new path for log
        profile_filename = timestamp + ".profile"
        return default_path / Path(profile_filename)

    # check if parent to a specified path exists
    if profile_path and not profile_path.parent.exists():
        logging.error("Path to output profile file dir is not valid")
        raise InvalidArgumentError("--profile_path", profile_path)

    if profile_path and profile_path.is_dir():
        logging.error("Path to output log file points to an existing folder")
        raise InvalidArgumentError("--log_path", profile_path)

    # otherwise we're good!
    return profile_path
