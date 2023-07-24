"""Module containing a class for a BiG-SCAPE input object, which has all the
input parameters/arguments"""

# from python
import logging
from enum import Enum
from pathlib import Path
from typing import Optional

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
        self.pfam_version: Optional[str] = None
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
        validate_pfam(self.pfam_version, self.pfam_path)
        validate_reference(self.mibig_version, self.reference_dir)


def validate_pfam(pfam_version, pfam_path):
    """Validates the pfam related properties"""

    # given only a path, the file must exist
    if not pfam_version and pfam_path and not pfam_path.exists():
        logging.error("Pfam file does not exist!")
        raise InvalidArgumentError("--pfam_path", pfam_path)

    # given a version (download action) the parent must exist
    if pfam_version and pfam_path and not pfam_path.parent.exists():
        logging.error("Path to pfam file is not valid")
        raise InvalidArgumentError("--pfam_path", pfam_path)

    if pfam_version and pfam_path and pfam_path.exists():
        logging.info(
            "Pfam file already exists, if you wish to re-download,"
            " delete or move this file. In the meantime, BiG_SCAPE"
            " will use the existing file."
        )


def validate_reference(mibig_version, reference_dir):
    """Validates the reference/MIBiG related properties"""

    # given neither -> do nothing

    # given reference dir and no mibig version, the dir must exist
    if not mibig_version and reference_dir and not reference_dir.exists():
        logging.error("GBK reference directory does not exist!")

    # given reference dir and mibig version (download action), the parent dir must exist
    if mibig_version and reference_dir and not reference_dir.parent.exists():
        logging.error("Path to pfam file is not valid")
        raise InvalidArgumentError("--reference_dir", reference_dir)

    if mibig_version and reference_dir and reference_dir.exists():
        logging.info(
            "GBK reference directory already exists, if you wish to "
            "re-download, delete or move this directory. In the "
            "meantime, BiG-SCAPE will use the existing directory."
        )


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
