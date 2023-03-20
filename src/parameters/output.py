"""Module containing a class for a BiG-SCAPE output object, which has all the
output parameters/arguments"""

# from python
from typing import Optional
from pathlib import Path

# from dependencies
# from other modules
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
        self.log_pat: Optional[Path] = None
        self.results: Optional[Path] = None
