"""Module containing a class for a BiG-SCAPE comparison object, which has all the
comparison parameters/arguments"""

# from other modules
from big_scape.errors import InvalidArgumentError

import big_scape.enums as bs_enums


class ComparisonParameters:
    """
    Class to store all run comparison parameters

    Attributes:
        alignment_mode: ALIGNMENT_MODE
        #TODO: proto/refion, similarity filters, distance calc modes, weights
    """

    def __init__(self) -> None:
        self.alignment_mode: bs_enums.ALIGNMENT_MODE = bs_enums.ALIGNMENT_MODE.AUTO

    def validate(self):
        """Validate the arguments contained in this object and set default values"""
        validate_alignment_mode(self.alignment_mode)


def validate_alignment_mode(alignment_mode) -> None:
    """Validate the passed alignment mode is one of the allowed modes"""
    # this should not ever run so long as the choices in the cmd_parser
    # remain updated and relevant
    valid_modes = [mode.value for mode in bs_enums.ALIGNMENT_MODE]
    if alignment_mode not in valid_modes:
        raise InvalidArgumentError("--align_mode", alignment_mode)
