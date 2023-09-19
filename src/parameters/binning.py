"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from other modules
import logging
from pathlib import Path
from typing import Optional
from src.errors import InvalidArgumentError


class BinningParameters:
    """Class to store all binning parameters

    Attributes:
        mix: bool
        legacy_bins: bool
        query_bgc_path: Optional[Path]
    """

    def __init__(self):
        self.mix: bool = False
        self.legacy_no_classify: bool = False
        self.query_bgc_path: Optional[Path] = None

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_mix(self.mix)
        validate_legacy_no_classify(self.legacy_no_classify)
        validate_query_bgc(self.query_bgc_path)

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

        # TODO: legacy_no_classify needs to be changed to legacy_classify, and add argument to no_classify
        if self.query_bgc_path is not None and (
            self.legacy_no_classify is False or self.mix is True
        ):
            logging.error(
                (
                    "The combination of arguments you have selected for binning is not valid."
                    " --query_bgc_path is not compatible with other binning arguments such"
                    " as --mix or --legacy_no_classify."
                )
            )
            raise InvalidArgumentError("--query_bgc_path", self.query_bgc_path)


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


def validate_query_bgc(query_bgc_path: Path):
    """Raises an InvalidArgumentError if the query bgc path does not exist"""

    if not query_bgc_path.exists():
        logging.error("Query BGC path does not exist!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)

    if not query_bgc_path.is_file():
        logging.error("Query BGC path is not a file!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)

    if query_bgc_path.suffix != ".gbk":
        logging.error("Query BGC path is not a .gbk file!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)
