"""Main file for BiG-SCAPE dereplicate module"""

# from python
import logging
from datetime import datetime


# from other modules
from big_scape.utility.version import get_bigscape_version
from big_scape.cli.config import BigscapeConfig
from big_scape.dereplicating.input_data_loading import load_input_data
from big_scape.dereplicating.sourmash_utilities import (
    make_sourmash_input,
    run_sourmash_branchwater,
    parse_sourmash_results,
)


def run_bigscape_dereplicate(run: dict) -> None:
    """Run bigscape dereplicate: clusters and dereplicates a set of BGCs"""

    bigscape_version = get_bigscape_version()

    logging.info(f"Running BiG-SCAPE {bigscape_version} dereplicate")

    logging.info("Using config file %s", run["config_file_path"])
    BigscapeConfig.parse_config(run["config_file_path"], run["log_path"])

    start_time = run["start_time"]

    # load input folder
    gbk_list = load_input_data(run)

    # generate sourmash input files
    sourmash_dir, cds_fasta_dir, manysketch_csv_path = make_sourmash_input(gbk_list, run)

    # run sourmash branchwater plugin in cmdline
    sourmash_pairwise_csv_path = run_sourmash_branchwater(run, sourmash_dir, cds_fasta_dir, manysketch_csv_path)

    # parse (sour)mash results
    edges = parse_sourmash_results(sourmash_pairwise_csv_path)

    # generate connected components & find cluster center
    # network = generate_network(edges)

    # write output
    # write_output(network, run)

    time_elapsed = datetime.now() - start_time

    logging.info(
        "BiG-SCAPE dereplicate finished in %s seconds", time_elapsed.total_seconds()
    )
