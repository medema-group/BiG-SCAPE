"""Module containing a class for a BiG-SCAPE hmmer object, which has all the
hmmer parameters/arguments"""

# from python
import logging
from pathlib import Path
from typing import Optional
import os

# from same module
from .run import RunParameters

# from other modules
from src.errors import InvalidArgumentError, ArgumentParseError


class HmmerParameters:
    """Class to store all HMMer parameters

    Attributes:
        force_hmmscan: bool
        force_hmmalign: bool
        domain_overlap_cutoff: float
        domain_includelist_path: Optional[Path]
    """

    def __init__(self):
        self.force_hmmscan: bool = False
        self.skip_hmmscan: bool = False
        self.domain_overlap_cutoff: float = 1.0
        self.domain_includelist_path: Optional[Path] = None

        # not from arguments:
        self.domain_includelist = []

    def validate(self, run: RunParameters):
        """Performs validation on this class parameters"""
        validate_hsp_overlap_cutoff(self.domain_overlap_cutoff)
        self.domain_includelist = validate_includelist(self.domain_includelist_path)
        validate_skip_hmmscan(self.skip_hmmscan, run.output.output_dir)


def validate_skip_hmmscan(skip_hmmscan: bool, output_dir: Path):
    """Validates whether an output directory exists and is not empty when running
    skip_hmm, which requires already processed gbk files and hance a DB in output"""

    if skip_hmmscan and not output_dir.exists():
        logging.error(
            "Output directory does not exist, skip_hmmscan requires a DB of\
                          already processed gbk files."
        )
        raise InvalidArgumentError("--skip_hmmscan", skip_hmmscan)

    if skip_hmmscan and output_dir.exists():
        contents = os.listdir(output_dir)
        if len(contents) == 0:
            logging.error(
                "Output directory empty, skip_hmmscan requires a DB of\
                          already processed gbk files."
            )
            raise InvalidArgumentError("--skip_hmmscan", skip_hmmscan)


def validate_includelist(
    domain_includelist_path: Optional[Path],
):
    """Validate the path to the domain include list and return a list of domain
    accession strings contained within this file

    Returns:
        list[str]: A list of domain accessions to include
    """

    # only validate if set
    if domain_includelist_path is None:
        return None

    if not domain_includelist_path.exists():
        logging.error("domain_includelist file does not exist!")
        raise InvalidArgumentError("--domain_includelist", domain_includelist_path)

    with domain_includelist_path.open(encoding="utf-8") as domain_includelist_file:
        lines = domain_includelist_file.readlines()

        lines = [line.strip() for line in lines]

        # expect Pfam accessions, i.e. PF00001 or PF00001.10
        lines_valid = map(
            lambda string: string.startswith("PF") and len(string) in range(7, 11),
            lines,
        )

        if not all(lines_valid):
            logging.error(
                "Invalid Pfam accession(s) found in file %s", domain_includelist_path
            )
            raise ArgumentParseError(
                "--domain_includelist_path", domain_includelist_path, ""
            )

        return lines


def validate_hsp_overlap_cutoff(cutoff: float):
    """Raises an InvalidArgumentError if cutoff is not between 0.0 and 1.0"""

    if cutoff < 0.0 or cutoff > 1.0:
        logging.error("Invalid cutoff (%f)! Must be between 0.0 and 1.0!", cutoff)
        raise InvalidArgumentError("--overlap_cutoff", cutoff)
