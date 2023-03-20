"""Module containing a class for a BiG-SCAPE hmmer object, which has all the
hmmer parameters/arguments"""

# from python
from typing import Optional, List

# from dependencies
# from other modules
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
