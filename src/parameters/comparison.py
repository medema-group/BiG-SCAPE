"""Module containing a class for a BiG-SCAPE comparison object, which has all the
comparison parameters/arguments"""

# from python
from typing import Optional


class ComparisonParameters:
    """
    Class to store all run comparison parameters

    Attributes:
        alignment_mode: str
        #TODO: proto/refion, similarity filters, distance calc modes, weigths
    """

    def __init__(self) -> None:
        self.alignment_mode: Optional[str] = None

    def parse(self, alignment_mode: str):
        """Load comparison arguments from commandline ArgParser object

        Args:
            alignment_mode (str): Alignment mode for each pair of gene clusters,
            choices=["global", "glocal", "auto"]
        """

        self.alignment_mode = alignment_mode