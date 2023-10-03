"""Contains errors raised for gene bank parsing"""


class InvalidGBKError(Exception):
    """Raised when a validation error occurs in parsing a GBK"""

    def __init__(self) -> None:
        super().__init__("Validation error occurred when parsing a GeneBank file")


class InvalidGBKRegionChildError(Exception):
    """Raised when a child to a GBK region is added that is not described in its parent
    qualifier
    """

    def __init__(self) -> None:
        super().__init__("GBK Feature is not parented correctly")
