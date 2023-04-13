"""Module containing a class for a BiG-SCAPE run networking object, which has all the
networking parameters/arguments"""

# from other modules
import logging
from src.errors import InvalidArgumentError


class NetworkingParameters:
    """
    Class to store all run networking parameters

    Attributes:
        gcf_cutoffs: List[float]
        include_singletons: Bool
        #TODO: modes
    """

    def __init__(self):
        self.gcf_cutoffs: list[float] = []
        self.include_singletons: bool = False

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        pass


def validate_gcf_cutoffs(gcf_cutoffs):
    """Raises an InvalidArgumentError if any cutoff is lower than 0"""
    for cutoff in gcf_cutoffs:
        if cutoff < 0:
            logging.error("One of the GCF cutoff values is invalid: %f < 0.00", cutoff)
            raise InvalidArgumentError("--gcf_cutoffs", gcf_cutoffs)
