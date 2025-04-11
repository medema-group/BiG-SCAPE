"""Main file for BiG-SCAPE dereplicate module"""

# from python
import logging
from datetime import datetime
import multiprocessing

# from other modules
from big_scape.utility.version import get_bigscape_version
from big_scape.cli.config import BigscapeConfig
from big_scape.dereplicating.input_data_loading import (
    load_input_folder,
    parse_gbk_files,
    gbk_factory,
)
import big_scape.enums as bs_enums
from big_scape.errors import InvalidGBKError


def run_bigscape_dereplicate(run: dict) -> None:
    """Run bigscape dereplicate: clusters and dereplicates a set of BGCs"""

    bigscape_version = get_bigscape_version()

    logging.info(f"Running BiG-SCAPE {bigscape_version} dereplicate")

    logging.info("Using config file %s", run["config_file_path"])
    BigscapeConfig.parse_config(run["config_file_path"], run["log_path"])

    start_time = run["start_time"]

    # load input folder

    input_gbk_files = load_input_folder(run)

    logging.info("Loading %d input GBKs", len(input_gbk_files))

    gbk_data = parse_gbk_files(input_gbk_files, bs_enums.SOURCE_TYPE.QUERY)

    # parse input GBKs

    cores = run["cores"]
    if cores is None:
        cores = multiprocessing.cpu_count()

    pool = multiprocessing.Pool(cores)

    data_package = map(lambda e: (e, run), gbk_data)

    gbk_list = pool.starmap(gbk_factory, data_package)

    if any(gbk is None for gbk in gbk_list):
        raise InvalidGBKError()

    logging.info("Succefully loaded %d input GBKs", len(gbk_list))

    # (log duplicated GBKs) and apply lenght constraints ???

    # concatenate CDSs

    # write CDS fastas

    # run (sour)mash

    # parse (sour)mash results

    # generate connected components

    # find cluster center

    # generate output

    time_elapsed = datetime.now() - start_time

    logging.info(
        "BiG-SCAPE dereplicate finished in %s seconds", time_elapsed.total_seconds()
    )
