"""Main benchmarking workflow"""

# from python
import sys

# from dependencies

# from other modules
from benchmark.parameters import parse_cmd
from benchmark.data import BenchmarkData
from benchmark.metrics import BenchmarkMetrics
from benchmark.output import OutputGenerator


def run_benchmark() -> None:
    """Run a benchmarking workflow. This will (if needed) run BiG-SCAPE on provided
    input data and evaluate computed GCFs based on a provided curated GCFs"""
    args = parse_cmd(sys.argv[1:])

    # load in both curated and copmuted GCF data
    db_path = args.computed_gcfs / "data_sqlite.db"
    data = BenchmarkData(args.curated_gcfs, db_path)
    data.load_curated_labels()
    data.load_computed_labels()

    # calculate metrics
    # v_measure = BenchmarkMetrics.calculate_v_measure(data)
    purities = BenchmarkMetrics.calculate_purity(data)
    entropies = BenchmarkMetrics.calculate_entropy(data)

    # output
    outputter = OutputGenerator(args.output_dir)
    outputter.initialize_output_dir()
    outputter.output_purities(purities)
    outputter.output_entropies(entropies)
