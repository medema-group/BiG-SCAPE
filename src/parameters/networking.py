"""Module containing a class for a BiG-SCAPE run networking object, which has all the
networking parameters/arguments"""

# from python
from typing import Optional

# from dependencies
# from other modules
# from this module


class Networking:
    """
    Class to store all run networking parameters

    Attributes:
        gcf_cutoff: float
        include_singletons: Bool
        #TODO: modes
    """

    def __init__(self) -> None:
        self.gcf_cutoff: Optional[float] = None
        self.include_ingletons: Optional[bool] = None
