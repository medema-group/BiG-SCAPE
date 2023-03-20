"""Module containing a class for a BiG-SCAPE comparison object, which has all the
comparison parameters/arguments"""

# from python
from typing import Optional

# from dependencies
# from other modules
# from this module


class Comparison:
    """
    Class to store all run comparison parameters

    Attributes:
        alignment_mode: str
        #TODO: proto/refion, similarity filters, distance calc modes, weigths
    """

    def __init__(self) -> None:
        self.alignment_mode: Optional[str] = None
