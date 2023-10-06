"""Main benchmarking workflow"""

# from python
import sys

# from dependencies

# from other modules
from benchmark.parameters import parse_cmd
from benchmark.data import BenchmarkData


def run_benchmark() -> None:
    """Run a benchmarking workflow. This will (if needed) run BiG-SCAPE on provided
    input data and evaluate computed GCFs based on a provided curated GCFs"""
    args = parse_cmd(sys.argv[1:])

    # load in both curated and copmuted GCF data
    data = BenchmarkData(args.curated_gcfs, args.computed_gcfs)
    data.load_curated_labels()
    data.load_computed_labels()

    # calculate metrics
