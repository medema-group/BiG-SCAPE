"""Module containing a class for a BiG-SCAPE diagnostics object, which has all the
diagnostics parameters/arguments"""

# from python
from typing import Optional

# from dependencies
# from other modules
# from this module


class Diagnostics:
    """
    Class to store all run diagnostics parameters

    Attributes:
        verbose: bool
        quiet: bool
        profiling: bool
        log_level: int
    """

    def __init__(self) -> None:
        self.verbose: Optional[bool] = None
        self.quiet: Optional[bool] = None
        self.profiling: Optional[bool] = None
        self.log_level: Optional[int] = None

    def parse(self, verbose: bool, quiet: bool, profiling: bool, log_level: int):
        """Load diagnostics arguments from commandline ArgParser object

        Args:
            verbose (bool): Output all kinds of logs
            quiet (bool): Do not print log info, just write to log file
            profiling (bool): Run profiler and output profile report
            log_level (int): Level of verbose intensity: 1,2,3
        """

        self.verbose = verbose
        self.quiet = quiet
        self.profiling = profiling
        self.log_level = log_level
