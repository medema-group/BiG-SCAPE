"""Main benchmarking workflow"""

# from python
import sys
import time

# from other modules
from benchmark.parameters import parse_cmd, validate_args
from benchmark.data import BenchmarkData
from benchmark.metrics import BenchmarkMetrics
from benchmark.output import OutputGenerator


def run_benchmark() -> None:
    """Benchmark: compare BiG-SCAPE output with curated GCF assignments"""
    args = parse_cmd(sys.argv[1:])
    validate_args(args)
    start_time = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    run_name = args.input_dir.stem + "_" + start_time

    # load in both curated and computed GCF data
    data = BenchmarkData(args.curated_gcfs, args.input_dir)
    data.load_curated_labels()
    data.load_computed_labels()

    # calculate metrics
    metrics = {}
    for fam_cutoff in data.computed_labels.keys():
        computed_labels_in_cutoff = data.computed_labels[fam_cutoff]
        calculator = BenchmarkMetrics(data.curated_labels, computed_labels_in_cutoff)
        metrics[fam_cutoff] = calculator.calculate_metrics()

        # output per cutoff
        metadata = OutputGenerator.generate_metadata(args, start_time, fam_cutoff)
        outputter = OutputGenerator(
            args.output_dir / f"cutoff_{fam_cutoff}", metadata, run_name
        )
        outputter.initialize_output_dir()
        outputter.output_purities(metrics[fam_cutoff]["purities"])
        outputter.output_entropies(metrics[fam_cutoff]["entropies"])
        outputter.output_matrix(metrics[fam_cutoff]["conf_matrix"])
        outputter.output_summary(
            metrics[fam_cutoff]["homogeneity"],
            metrics[fam_cutoff]["completeness"],
            metrics[fam_cutoff]["v_measure"],
            metrics[fam_cutoff]["purities"],
            metrics[fam_cutoff]["entropies"],
            metrics[fam_cutoff]["associations"],
            metrics[fam_cutoff]["summary_stats"],
        )
        outputter.plot_conf_matrix_heatmap(metrics[fam_cutoff]["conf_matrix"])
    # output summary per cutoff
    metadata = OutputGenerator.generate_metadata(args, start_time)
    outputter = OutputGenerator(args.output_dir, metadata, run_name)
    outputter.plot_per_cutoff(metrics)
    outputter.output_summary_per_cutoff(metrics)
