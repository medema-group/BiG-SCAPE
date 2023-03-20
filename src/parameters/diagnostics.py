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
