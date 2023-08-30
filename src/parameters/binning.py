"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from other modules
import logging
from src.errors import InvalidArgumentError


class BinningParameters:
    """Class to store all binning parameters

    Attributes:
        mix: bool
        legacy_bins: bool
    """

    def __init__(self):
        self.mix: bool = False
        self.legacy_no_classify: bool = False

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_mix(self.mix)
        validate_legacy_no_classify(self.legacy_no_classify)

        self.validate_work()

    def validate_work(self):
        """Raise an error if the combination of parameters in this object means no work will be done"""

        if self.mix is False and self.legacy_no_classify is True:
            logging.error(
                (
                    "The combination of arguments you have selected for binning means no work will "
                    "be done. Please add either --mix, or remove --legacy_no_classify in order to enable "
                    "comparisons"
                )
            )
            raise InvalidArgumentError("--mix", self.mix)


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


def validate_legacy_no_classify(legacy_no_classify: bool):
    if legacy_no_classify is None:
        logging.error(
            (
                "'--legacy_no_classify' is somehow set to None. If you did not make any changes to "
                "the code, contact the maintainers by submitting your dataset and an "
                "issue."
            )
        )
        raise InvalidArgumentError("--legacy_no_classify", legacy_no_classify)
