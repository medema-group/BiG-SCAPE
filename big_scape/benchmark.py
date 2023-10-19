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
    for fam_cutoff in data.computed_labels.keys():
        computed_labels_in_cutoff = data.computed_labels[fam_cutoff]
        calculator = BenchmarkMetrics(data.curated_labels, computed_labels_in_cutoff)
        metrics[fam_cutoff] = calculator.calculate_metrics()

        # output per cutoff
        logging.info("Generating cutoff %s output", fam_cutoff)
        metadata = OutputGenerator.generate_metadata(run, fam_cutoff)
        outputter = OutputGenerator(
            run["output_dir"] / f"cutoff_{fam_cutoff}", metadata, run["label"]
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
    logging.info("Generating summary output")
    # output summary per cutoff
    metadata = OutputGenerator.generate_metadata(run)
    outputter = OutputGenerator(run["output_dir"], metadata, run["label"])
    outputter.plot_per_cutoff(metrics)
    outputter.output_summary_per_cutoff(metrics)
    logging.info("Benchmark done!")
    # print("running bigscape benchmark")
    # print(run)
