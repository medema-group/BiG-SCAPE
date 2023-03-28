"""Module containing a class for a BiG-SCAPE run binning object, which has all the
binning parameters/arguments"""

# from python
from typing import Optional


class BinningParameters:
    """
    Class to store all run binning parameters

    Attributes:
        mix: Bool
        #TODO: modes, banned_classes, no_classify
    """

    def __init__(self) -> None:
        self.mix: Optional[bool] = None

    def parse(self, mix: bool):
        """Load networking arguments from commandline ArgParser object

        Args:
            mix (bool): Run an all-vs-all analysis
        """

        self.mix = mix
