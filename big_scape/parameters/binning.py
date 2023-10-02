"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from python
import logging
from pathlib import Path
from typing import Optional

# from other modules
from big_scape.errors import InvalidArgumentError


class BinningParameters:
    """Class to store all binning parameters

    Attributes:
        mix: bool
        legacy_bins: bool
        query_bgc_path: Optional[Path]
    """

    def __init__(self):
        self.no_mix: bool = False
        self.legacy_classify: bool = False
        self.classify: bool = False
        self.query_bgc_path: Optional[Path] = None

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_mix(self.no_mix)
        validate_legacy_classify(self.legacy_classify)
        validate_query_bgc(self.query_bgc_path)

        self.validate_work()

    def validate_work(self):
        """Raise an error if the combination of parameters in this object means no work
        will be done
        """

        if (
            self.no_mix is True
            and self.legacy_classify is False
            and self.classify is False
        ):
            logging.error(
                (
                    "The combination of arguments you have selected for binning means no work will "
                    "be done. Please remove either --no_mix, or add --legacy_no_classify/--classify"
                    " in order to enable comparisons"
                )
            )
            raise InvalidArgumentError("--no_mix", self.no_mix)

        # TODO: legacy_no_classify needs to be changed to legacy_classify, and add
        # argument to no_classify
        if self.query_bgc_path is not None and (
            self.legacy_classify is True or self.no_mix is True or self.classify is True
        ):
            logging.error(
                (
                    "The combination of arguments you have selected for binning is not valid."
                    " --query_bgc_path is not compatible with other binning arguments such"
                    " as --no_mix, --classify or --legacy_classify."
                )
            )
            raise InvalidArgumentError("--query_bgc_path", self.query_bgc_path)


def validate_mix(no_mix: bool):
    """Validates the no_mix attribute. Raises an exception if it is set to none, which
    should never happen
    """
    if no_mix is None:
        logging.error(
            (
                "'--mix' is somehow set to None. If you did not make any changes to "
                "the code, contact the maintainers by submitting your dataset and an "
                "issue."
            )
        )
        raise InvalidArgumentError("--no_mix", no_mix)


def validate_legacy_classify(legacy_classify: bool):
    """Validates the legacy_classify attribute. Raises an exception if it is set to
    none, which should never happen
    """
    if legacy_classify is None:
        logging.error(
            (
                "'--legacy_classify' is somehow set to None. If you did not make any changes to "
                "the code, contact the maintainers by submitting your dataset and an "
                "issue."
            )
        )
        raise InvalidArgumentError("--legacy_classify", legacy_classify)


def validate_classify(classify: bool):
    """Validates the classify attribute. Raises an exception if it is set to none, which
    should never happen
    """
    if classify is None:
        logging.error(
            (
                "'--classify' is somehow set to None. If you did not make any changes to "
                "the code, contact the maintainers by submitting your dataset and an "
                "issue."
            )
        )
        raise InvalidArgumentError("--classify", classify)


def validate_query_bgc(query_bgc_path: Path):
    """Raises an InvalidArgumentError if the query bgc path does not exist"""

    if query_bgc_path and not query_bgc_path.exists():
        logging.error("Query BGC path does not exist!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)

    if query_bgc_path and not query_bgc_path.is_file():
        logging.error("Query BGC path is not a file!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)

    if query_bgc_path and query_bgc_path.suffix != ".gbk":
        logging.error("Query BGC path is not a .gbk file!")
        raise InvalidArgumentError("--query_bgc_path", query_bgc_path)
