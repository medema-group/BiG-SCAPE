"""Contains the run object, which details the structure of the namespace as parsed in
cmd_parser
"""


# from python
from datetime import datetime
from argparse import Namespace
from multiprocessing import cpu_count

# from other modules
from src.diagnostics import init_logger, init_logger_file

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
        legacy: bool
        input: src.parameters.InputParameters
        hmmer: src.parameters.HmmerParameters
        binning: src.parameters.BinningParameters
        comparison: src.parameters.ComparisonParameters
        networking: src.parameters.NetworkingParameters
        diagnostics: src.parameters.DiagnosticsParameters
        output: src.parameters.OutputParameters
    """

    def __init__(self):
        self.label: str = "BiG-SCAPE"
        self.cores: int = cpu_count()
        self.legacy: bool = False
        self.input = InputParameters()
        self.hmmer = HmmerParameters()
        self.binning = BinningParameters()
        self.comparison = ComparisonParameters()
        self.networking = NetworkingParameters()
        self.diagnostics = DiagnosticsParameters()
        self.output = OutputParameters()

    def start(self) -> datetime:
        """Start this run, set the label and return a datetime of the start time

        Returns:
            datetime: datetime object corresponding to the start of the run
        """
        start_time: datetime = datetime.now()

        timestamp = start_time.strftime("%d-%m-%Y %H_%M_%S")
        self.label = f"{self.label}_{timestamp}"

        return start_time

    def validate(self):
        """Executes validation on everything, setting default values and returning
        errors if anything is wrong
        Also initializes the logger
        """
        # diagnostics validate must not use logger
        self.diagnostics.validate()

        # initializes the logger, needed for the validations below
        init_logger(self)

        self.input.validate()
        self.hmmer.validate(self)
        self.binning.validate()
        self.comparison.validate()
        self.networking.validate()
        self.output.validate()
        self.hmmer.validate(self)
        self.binning.validate()
        self.comparison.validate()
        self.diagnostics.validate()
        self.networking.validate()

        # initializes the logger file, can only happen once output is validated
        # and the log file path is created
        init_logger_file(self)

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
