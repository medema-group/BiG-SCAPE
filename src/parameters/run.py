"""Contains the run object, which details the structure of the namespace as parsed in
cmd_parser
"""


# from python
from argparse import Namespace
from multiprocessing import cpu_count

# from this module
from src.parameters.input import InputParameters
from src.parameters.hmmer import HmmerParameters
from src.parameters.binning import BinningParameters
from src.parameters.comparison import ComparisonParameters
from src.parameters.networking import NetworkingParameters
from src.parameters.diagnostics import DiagnosticsParameters
from src.parameters.output import OutputParameters


class RunParameters(Namespace):
    """
    Class to store all run parameters

    Attributes:
        label: str
        cores: int
        input: Input, object storing all
        output: Output
        diagnostics: Diagnostics
        hmmer: Hmmer
        networking: Networking
        binning: Binning
        comparison: Comparison
    """

    def __init__(self):
        self.label: str = ""
        self.cores: int = cpu_count()
        self.input = InputParameters()
        self.hmmer = HmmerParameters()
        self.binning = BinningParameters()
        self.comparison = ComparisonParameters()
        self.networking = NetworkingParameters()
        self.diagnostics = DiagnosticsParameters()
        self.output = OutputParameters()

    def validate(self):
        """Executes validation on everything, setting default values and returning
        errors if anything is wrong
        """
        self.input.validate()

        self.output.validate()

    # from https://stackoverflow.com/questions/18668227
    # this method is called when a new attribute is set on this class. specifically,
    # this is called when argparser tries to add an attribute like "input.gbk_dir",
    # including the dot. this method ensures that instead a property will be added on
    # the input class instead
    def __setattr__(self, name, value):
        if "." in name:
            group, name = name.split(".", 1)
            ns = getattr(self, group, RunParameters())
            # this part added to check if property exists on object
            if name not in ns.__dict__:
                raise KeyError(
                    (
                        f"parameter subgroup '{group}' does not contain property "
                        f"'{name}'! "
                        "This usually means one of the 'dest' arguments to "
                        "parser.add_argument is wrong or the parameter is not present "
                        "in a parameter subclass"
                    )
                )
            setattr(ns, name, value)
            self.__dict__[group] = ns
        else:
            self.__dict__[name] = value
