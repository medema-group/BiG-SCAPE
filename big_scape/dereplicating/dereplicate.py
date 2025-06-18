"""Main file for BiG-SCAPE dereplicate module"""

# from python
import logging
from datetime import datetime

# from other modules
from big_scape.utility.version import get_bigscape_version
from big_scape.dereplicating.input_data_loading import load_input_data
from big_scape.dereplicating.sourmash_utilities import (
    make_sourmash_input,
    run_sourmash_branchwater,
    parse_sourmash_results,
)
from big_scape.dereplicating.output_generation import write_output
from big_scape.dereplicating.networking import Network


def run_bigscape_dereplicate(run: dict) -> None:
    """Run bigscape dereplicate: clusters and dereplicates a set of BGCs"""

    bigscape_version = get_bigscape_version()

    logging.info(f"Running BiG-SCAPE {bigscape_version} dereplicate")

    start_time = run["start_time"]

    # load input folder
    logging.info("Loading input data")
    gbk_list = load_input_data(run)

    # generate sourmash input files
    logging.info("Generating sourmash input files")
    sourmash_dir, cds_fasta_dir, manysketch_csv_path = make_sourmash_input(gbk_list, run)

    # run sourmash branchwater plugin in cmdline
    logging.info("Running sourmash (with the branchwater plugin)")
    sourmash_pairwise_csv_path = run_sourmash_branchwater(run, sourmash_dir, cds_fasta_dir, manysketch_csv_path)

    # parse (sour)mash results
    logging.info("Parsing sourmash results")
    edges, nodes = parse_sourmash_results(sourmash_pairwise_csv_path, run['cutoff'])

    # generate connected components & find cluster center
    logging.info("Generating connected components and finding representative clusters")
    network = Network(edges, nodes)

    # write output
    logging.info("Writing output files")
    write_output(network, run, gbk_list)

    time_elapsed = datetime.now() - start_time

    logging.info(
        "BiG-SCAPE dereplicate finished in %s seconds", time_elapsed.total_seconds()
    )
