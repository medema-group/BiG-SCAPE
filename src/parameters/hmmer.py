"""Module containing a class for a BiG-SCAPE hmmer object, which has all the
hmmer parameters/arguments"""

# from python
import logging
from typing import Optional, List
from pathlib import Path

# from dependencies
# from other modules
from src.errors import InvalidInputArgError

# from this module


class Hmmer:
    """
    Class to store all run hmmer parameters

    Attributes:
        domain_overlap_cutoff: float
        force_hmmscan: Bool
        skip_alignment: Bool
        domain_includelist: List[str]
    """

    def __init__(self) -> None:
        self.domain_overlap_cutoff: Optional[float] = None
        self.force_hmmscan: Optional[bool] = None
        self.skip_aligmnent: Optional[bool] = None
        self.domain_includelist: Optional[List[str]] = None

    def parse(
        self,
        domain_overlap_cutoff: float,
        force_hmmscan: bool,
        skip_alignment: bool,
        domain_includelist_path: Path,
    ):
        """Load hmmer arguments from commandline ArgParser object

        Args:
            domain_overlap_cutoff (float): overlap percentage at which domains are
            considered to overlap
            force_hmmscan (bool): Force domain prediction
            skip_alignment (bool): Skip multiple alignment
            domain_includelist_path (Path): Path to txt file with Pfam accessions
        """

        self.domain_overlap_cutoff = domain_overlap_cutoff
        self.force_hmmscan = force_hmmscan
        self.skip_aligmnent = skip_alignment

        if domain_includelist_path and not domain_includelist_path.exists():
            logging.error("Path to domain_includelist file is not valid")
            raise InvalidInputArgError()

        elif not domain_includelist_path:
            pass

        else:
            with domain_includelist_path.open(
                encoding="utf-8"
            ) as domain_includelist_file:
                lines = domain_includelist_file.readlines()

                lines = [line.strip() for line in lines]

                # expect Pfam accessions, i.e. PF00001 or PF00001.10
                lines_valid = map(
                    lambda string: string.startswith("PF")
                    and len(string) in range(7, 11),
                    lines,
                )

                if not all(lines_valid):
                    logging.error("Invalid Pfam accession(s)")
                    raise InvalidInputArgError

                self.domain_includelist = lines

                # TODO: test this
