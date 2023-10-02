"""Module containing a class for a BiG-SCAPE comparison object, which has all the
comparison parameters/arguments"""

# from python
from enum import Enum

# from other modules
from big_scape.errors import InvalidArgumentError


class ALIGNMENT_MODE(Enum):
    GLOBAL = "global"
    GLOCAL = "glocal"
    AUTO = "auto"


class ComparisonParameters:
    """
    Class to store all run comparison parameters

    Attributes:
        alignment_mode: ALIGNMENT_MODE
        #TODO: proto/refion, similarity filters, distance calc modes, weights
    """

    def __init__(self):
        self.alignment_mode: ALIGNMENT_MODE = ALIGNMENT_MODE.AUTO

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_alignment_mode(self.alignment_mode)


def validate_alignment_mode(alignment_mode):
    """Validate the passed alignment mode is one of the allowed modes"""
    # this should not ever run so long as the choices in the cmd_parser
    # remain updated and relevant
    valid_modes = [mode.value for mode in ALIGNMENT_MODE]
    if alignment_mode not in valid_modes:
        raise InvalidArgumentError("--align_mode", alignment_mode)
