"""Contains errors raised for input argument parsing"""

# from python
from typing import Any


class InvalidArgumentError(Exception):
    """Raised when a validation error occurs in parsing a commandlind input
    argument"""

    def __init__(self, argument_name: str, argument_value: Any):
        super().__init__(
            (
                f"Validation error occured for argument {argument_name} "
                f"with value: {argument_value}. ",
                "See bigscape.py --help for usage information",
            )
        )


class ArgumentParseError(Exception):
    """Raised when a parameter is valid, but after parsing/validating the content
    pointed to by the parameter, some error occurs
    """

    def __init__(self, argument_name, argument_value, description) -> None:
        super().__init__(
            (
                f"Parsing error occurred for argument {argument_name} ",
                f"with value: {argument_value}. {description}",
            )
        )
