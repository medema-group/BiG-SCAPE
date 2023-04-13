"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from other modules
import logging
from src.errors import InvalidArgumentError


class BinningParameters:
    """Class to store all binning parameters

    Attributes:
        mix: bool
    """

    def __init__(self):
        self.mix: bool = True

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_mix(self.mix)


def validate_mix(mix: bool):
    """Validates the mix attribute. Raises an exception if it is set to none, which
    should never happen
    """
    if mix is None:
        logging.error(
            (
                "'--mix' is somehow set to None. If you did not make any changes to "
                "the code, contact the maintainers by submitting your dataset and an "
                "issue."
            )
        )
        raise InvalidArgumentError("--mix", mix)
