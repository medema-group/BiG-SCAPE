"""Module containing a class for a BiG-SCAPE run object, which has all the
run parameters/arguments"""

# from python
from typing import Optional

# from dependencies
# from other modules
# from this module
from src.parameters.input import Input
from src.parameters.output import Output
from src.parameters.comparison import Comparison
from src.parameters.diagnostics import Diagnostics
from src.parameters.hmmer import Hmmer
from src.parameters.binning import Binning
from src.parameters.networking import Networking


class Run:
    """
    Class to store all run parameters

    Attributes:
        input: Input, object storing all
        output: Output
        label: str
        cores: int
        diagnostics: Diagnostics
        hmmer: Hmmer
        networking: Networking
        binning: Binning
        comparison: Comparison
    """

    def __init__(self) -> None:
        self.input: Optional[Input] = None
        self.output: Optional[Output] = None
        self.label: Optional[str] = None
        self.cores: Optional[int] = None
        self.diagnostics: Optional[Diagnostics] = None
        self.hmmer: Optional[Hmmer] = None
        self.networking: Optional[Networking] = None
        self.binning: Optional[Binning] = None
        self.comparison: Optional[Comparison] = None
