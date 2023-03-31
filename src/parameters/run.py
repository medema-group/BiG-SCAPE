"""Module containing a class for a BiG-SCAPE run object, which has all the
run parameters/arguments"""

# from python
from typing import Optional, Any

# from this module
from src.parameters.input import InputParameters
from src.parameters.output import OutputParameters
from src.parameters.comparison import ComparisonParameters
from src.parameters.diagnostics import DiagnosticsParameters
from src.parameters.hmmer import HmmerParameters
from src.parameters.binning import BinningParameters
from src.parameters.networking import NetworkingParameters


class RunParameters:
    """
    Class to store all run parameters

    Attributes:
        input: Input, object storing all
        output: Output
        label: str
        cores: int
        diagnostics: Diagnostics
        hmmer: Hmmer
        networking: Networking
        binning: Binning
        comparison: Comparison
    """

    def __init__(self) -> None:
        self.input: Optional[InputParameters] = InputParameters()
        self.output: Optional[OutputParameters] = OutputParameters()
        self.label: Optional[str] = None
        self.cores: Optional[int] = None
        self.diagnostics: Optional[DiagnosticsParameters] = DiagnosticsParameters()
        self.hmmer: Optional[HmmerParameters] = HmmerParameters()
        self.networking: Optional[NetworkingParameters] = NetworkingParameters()
        self.binning: Optional[BinningParameters] = BinningParameters()
        self.comparison: Optional[ComparisonParameters] = ComparisonParameters()

    def parse(self, args: Any):
        """Load input arguments from commandline ArgParser object

        Args:
            args (any): object storing input arguments
        """

        self.cores = args.cores
        self.label = args.label

        if self.input is None:
            raise ValueError()
        self.input.parse(
            args.inputdir,
            args.input_mode,
            args.dataset_path,
            args.metadata_path,
            args.pfam_path,
            args.download_pfam,
            args.mibig_version,
            args.reference_dir,
            args.include_gbk,
            args.exclude_gbk,
            args.query_bgc_path,
            args.min_bgc_length,
        )

        if self.output is None:
            raise ValueError()
        self.output.parse(args.db_path, args.log_path, args.outputdir)

        if self.diagnostics is None:
            raise ValueError()
        self.diagnostics.parse(args.verbose, args.quiet, args.profiling)

        if self.hmmer is None:
            raise ValueError()
        self.hmmer.parse(
            args.domain_overlap_cutoff,
            args.force_hmmscan,
            args.force_hmmalign,
            args.domain_includelist_path,
        )

        if self.networking is None:
            raise ValueError()
        self.networking.parse(args.gcf_cutoffs, args.include_singletons)

        if self.binning is None:
            raise ValueError()
        self.binning.parse(args.mix)

        if self.comparison is None:
            raise ValueError()
        self.comparison.parse(args.alignment_mode)
