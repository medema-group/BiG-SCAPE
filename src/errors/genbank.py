"""Contains errors raised for gene bank parsing"""


class InvalidGBKError(Exception):
    """Raised when a validation error occurs in parsing a GBK"""

    def __init__(self):
        super().__init__("Validation error occurred when parsing a GeneBank file")
