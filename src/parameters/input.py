"""Module containing a class for a BiG-SCAPE input object, which has all the
input parameters/arguments"""

# from python
import logging
from enum import Enum
from pathlib import Path
from typing import Optional
import os

# from other modules
from src.errors import InvalidArgumentError


class INPUT_MODE(Enum):
    FLAT = "flat"
    RECURSIVE = "recursive"


class InputParameters:
    """Class to store all data input related parameters

    Attributes:
        input_dir: Path
        input_mode: INPUT_MODE
        pfam_path: Optional[Path]
        pfam_version: bool
        mibig_dir: Optional[Path]
        mibig_version: str
        metadata_path: Optional[Path]
        dataset_path: Optional[Path]
        reference_dir: Optional[Path]
        query_bgc_path: Optional[Path]
        include_gbk: list[str]
        exclude_gbk: list[str]
        min_bgc_length: int
    """

    def __init__(self):
        # main input
        self.input_dir: Path = Path("")
        self.input_mode: INPUT_MODE = INPUT_MODE.RECURSIVE

        # other db
        # TODO: test
        self.pfam_path: Optional[Path] = None
        self.mibig_version: Optional[str] = None

        # other user supplied stuff
        # TODO: test
        self.metadata_path: Optional[Path] = None
        self.dataset_path: Optional[Path] = None
        self.reference_dir: Optional[Path] = None
        self.query_bgc_path: Optional[Path] = None

        # filtering
        # TODO: test
        self.include_gbk: list[str] = []
        self.exclude_gbk: list[str] = []
        self.min_bgc_length: int = 0
        # TODO: this one has a test
        self.cds_overlap_cutoff: Optional[float] = None

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_input_dir(self.input_dir)
        validate_input_mode(self.input_mode)
        validate_pfam(self.pfam_path)
        validate_reference(self.reference_dir)
        validate_query_bgc(self.query_bgc_path)
        validate_cds_overlap_cutoff(self.cds_overlap_cutoff)


def validate_pfam(pfam_path):
    """Validates the pfam related properties"""

    # given only a path, the file must exist
    if pfam_path and not pfam_path.exists():
        logging.error("Pfam file does not exist!")
        raise InvalidArgumentError("--pfam_path", pfam_path)


def validate_reference(reference_dir):
    """Validates the reference/MIBiG related properties"""

    # given reference dir, the dir must exist
    if reference_dir and not reference_dir.exists():
        logging.error("GBK reference directory does not exist!")
        raise InvalidArgumentError("--reference_dir", reference_dir)

    if reference_dir and reference_dir.exists():
        contents = os.listdir(reference_dir)
        if len(contents) == 0:
            logging.error("GBK reference directory empty!")
            raise InvalidArgumentError("--reference_dir", reference_dir)


def validate_query_bgc(query_bgc_path):
    """Validates the query_gbc_path property"""

    # given single query bgc file, file must exist
    if query_bgc_path and not query_bgc_path.exists():
        logging.error("Query BGC file does not exist!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)

    if query_bgc_path and not query_bgc_path.is_file():
        logging.error("Query BGC file is not a file!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)


def validate_input_dir(input_dir):
    """Validates the gbk_dir property"""
    if input_dir is None:
        logging.error("GBK Input directory is not set!")
        raise InvalidArgumentError("--input_dir", input_dir)

    if not input_dir.exists():
        logging.error("GBK Input directory does not exist!")
        raise InvalidArgumentError("--input_dir", input_dir)

    if not input_dir.is_dir():
        logging.error("GBK Input directory is not a directory!")
        raise InvalidArgumentError("--input_dir", input_dir)


def validate_input_mode(input_mode):
    """validates the input_mode property. Raises an InvalidArgumentError if the
    input_mode parameter is invalid
    """
    # check if the property matches one of the enum values
    valid_modes = [mode.value for mode in INPUT_MODE]
    matches = input_mode in valid_modes
    if not matches:
        logging.error("Invalid input mode. Must be of type: %s", ", ".join(valid_modes))
        raise InvalidArgumentError("--input_mode", input_mode)


def validate_cds_overlap_cutoff(cutoff: float):
    """Raises an InvalidArgumentError if cutoff is not between 0.0 and 1.0"""

    if cutoff < 0.0 or cutoff > 1.0:
        logging.error("Invalid cutoff (%f)! Must be between 0.0 and 1.0!", cutoff)
        raise InvalidArgumentError("--overlap_cutoff", cutoff)
