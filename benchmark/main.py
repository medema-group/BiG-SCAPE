"""Main benchmarking workflow"""

# from python
import sys

# from other modules
from benchmark.parameters import parse_cmd, validate_args
from benchmark.data import BenchmarkData
from benchmark.metrics import BenchmarkMetrics
from benchmark.output import OutputGenerator


def run_benchmark() -> None:
    """Benchmark: compare BiG-SCAPE output with curated GCF assignments"""
    args = parse_cmd(sys.argv[1:])
    validate_args(args)

    # load in both curated and computed GCF data
    db_path = args.bigscape_dir / "data_sqlite.db"
    data = BenchmarkData(args.curated_gcfs, db_path)
    data.load_curated_labels()
    data.load_computed_labels()

    # calculate metrics
    for fam_cutoff in data.computed_labels.keys():
        computed_labels_in_cutoff = data.computed_labels[fam_cutoff]
        metrics = BenchmarkMetrics(data.curated_labels, computed_labels_in_cutoff)
        homogeneity, completeness, v_measure = metrics.calculate_v_measure()
        purities = metrics.calculate_purity()
        entropies = metrics.calculate_entropy()
        associations = metrics.compare_association()
        summary_stats = metrics.calculate_summary()

        # output
        outputter = OutputGenerator(args.output_dir / f"cutoff_{fam_cutoff}")
        outputter.initialize_output_dir()
        outputter.output_purities(purities)
        outputter.output_entropies(entropies)
        outputter.output_summary(
            homogeneity,
            completeness,
            v_measure,
            purities,
            entropies,
            associations,
            summary_stats,
        )
