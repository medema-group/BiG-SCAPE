"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from python
from typing import Optional

# from dependencies
# from other modules
# from this module


class Binning:
    """
    Class to store all run binning parameters

    Attributes:
        mix: Bool
        #TODO: modes, banned_classes, no_classify
    """

    def __init__(self) -> None:
        self.mix: Optional[bool] = None
