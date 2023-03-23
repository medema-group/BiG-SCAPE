"""Contains errors raised for input argument parsing"""


class InvalidInputArgError(Exception):
    """Raised when a validation error occurs in parsing a commandlind input
    argument"""

    def __init__(self):
        super().__init__(
            "Validation error occurred when parsing a command\
                         line input argument"
        )
