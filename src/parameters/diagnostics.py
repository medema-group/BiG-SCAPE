"""Module containing a class for a BiG-SCAPE diagnostics object, which has all the
diagnostics parameters/arguments"""

# from python
from typing import Optional


class DiagnosticsParameters:
    """
    Class to store all run diagnostics parameters

    Attributes:
        verbose: bool
        quiet: bool
        profiling: bool
    """

    def __init__(self) -> None:
        self.verbose: Optional[bool] = None
        self.quiet: Optional[bool] = None
        self.profiling: Optional[bool] = None

    def parse(self, verbose: bool, quiet: bool, profiling: bool):
        """Load diagnostics arguments from commandline ArgParser object

        Args:
            verbose (bool): Output all kinds of logs
            quiet (bool): Do not print log info, just write to log file
        """

        self.verbose = verbose
        self.quiet = quiet
        self.profiling = profiling
