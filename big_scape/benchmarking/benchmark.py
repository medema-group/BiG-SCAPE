"""Main file for BiG-SCAPE benchmark module"""

# from python
from typing import Any
import logging

# from other modules
from big_scape.benchmarking import BenchmarkData
from big_scape.benchmarking import BenchmarkMetrics
from big_scape.benchmarking import OutputGenerator


def run_bigscape_benchmark(run: dict):
    """Run bigscape benchmark: compares the results of a query against
    a benchmark set of BGC <-> GCF assignments"""

    data = BenchmarkData(run["gcf_assignment_file"], run["big_dir"])
    data.load_curated_labels()
    data.load_computed_labels()

    logging.info("Calculating benchmarking metrics per cutoff")
    metrics: dict[str, dict[str, Any]] = {}
    for fam_cutoff in sorted(data.computed_labels.keys(), key=float):
        computed_labels_in_cutoff = data.computed_labels[fam_cutoff]
        calculator = BenchmarkMetrics(data.curated_labels, computed_labels_in_cutoff)
        metrics[fam_cutoff] = calculator.calculate_metrics()

        # output per cutoff
        logging.info("Generating cutoff %s output", fam_cutoff)
        metadata = OutputGenerator.generate_metadata(run, fam_cutoff)
        cutoff_path = run["output_dir"] / f"cutoff_{fam_cutoff}"
        outputter = OutputGenerator(cutoff_path, metadata, run["label"])
        outputter.output_metrics(metrics, fam_cutoff)

    logging.info("Generating summary output")
    # output summary of all cutoffs
    metadata = OutputGenerator.generate_metadata(run)
    outputter = OutputGenerator(run["output_dir"], metadata, run["label"], data.tool)
    outputter.plot_per_cutoff(metrics)
    outputter.output_summary_per_cutoff(metrics)
    logging.info("Benchmark done!")
