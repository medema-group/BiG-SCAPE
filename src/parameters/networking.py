"""Module containing a class for a BiG-SCAPE run networking object, which has all the
networking parameters/arguments"""

# from python
import logging
from typing import Optional, List

# from other modules
from src.errors import InvalidInputArgError


class NetworkingParameters:
    """
    Class to store all run networking parameters

    Attributes:
        gcf_cutoffs: List[float]
        include_singletons: Bool
        #TODO: modes
    """

    def __init__(self) -> None:
        self.gcf_cutoffs: Optional[List[float]] = None
        self.include_ingletons: Optional[bool] = None

    def parse(self, gcf_cutoffs: str, include_singletons: bool):
        """Load networking arguments from commandline ArgParser object

        Args:
            gcf_cutoffs (string): Distance cutoff values
            include_ingletons (bool): Include nodes that have no edges
        """

        self.include_ingletons = include_singletons

        gcf_cutoffs_list = gcf_cutoffs.split(",")
        try:
            gcf_cutoffs_list_float = [float(cutoff) for cutoff in gcf_cutoffs_list]
            self.gcf_cutoffs = gcf_cutoffs_list_float
        except ValueError:
            logging.error("value given for gcf_cutoffs invalid")
            raise InvalidInputArgError()
