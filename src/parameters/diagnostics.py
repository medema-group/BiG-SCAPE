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

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        pass
