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
        v_measure = BenchmarkMetrics.calculate_v_measure(
            data.curated_labels, computed_labels_in_cutoff
        )
        purities = BenchmarkMetrics.calculate_purity(
            data.curated_labels, computed_labels_in_cutoff
        )
        entropies = BenchmarkMetrics.calculate_entropy(
            data.curated_labels, computed_labels_in_cutoff
        )
        (
            nr_cur_fams,
            nr_comp_fams,
            avg_cur_size,
            avg_comp_size,
        ) = BenchmarkMetrics.calculate_summary(
            data.curated_labels, computed_labels_in_cutoff
        )

        # output
        outputter = OutputGenerator(args.output_dir / f"cutoff_{fam_cutoff}")
        outputter.initialize_output_dir()
        outputter.output_v_measure(v_measure)
        outputter.output_purities(purities)
        outputter.output_entropies(entropies)
        outputter.output_summary(
            v_measure,
            purities,
            entropies,
            nr_cur_fams,
            nr_comp_fams,
            avg_cur_size,
            avg_comp_size,
        )
