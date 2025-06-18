""" Click parameters for the BiG-SCAPE Dereplicate CLI command """

# from python
import click
import logging

# from other modules
from big_scape.dereplicating.dereplicate import run_bigscape_dereplicate
from big_scape.diagnostics import init_logger, init_logger_file
import big_scape.enums as bs_enums
from big_scape.cli.config import BigscapeConfig

# from this module
from .cli_common_options import common_all, common_cluster_query_dereplicate
from .cli_validations import (
    set_start,
    validate_output_paths,
)


@click.command()
@common_all
@common_cluster_query_dereplicate
@click.option(
    "--cutoff",
    type=click.FloatRange(min=0.0, max=1.0),
    default=0.8,
    help=(
        "Similarity threshold for sourmash distances. Only pairs with a "
        "similarity equal or above this value will be considered for clustering."
    ),
)
@click.pass_context
def dereplicate(ctx, *args, **kwargs):
    """
    BiG-SCAPE - DEREPLICATE

    Dereplicate mode - BiG-SCAPE performs a pairwise comparison
    of BGCs based on the protein sequence comparison tool sourmash,
    clusters them based on a similarity threshold, and selects
    a representative BGC per cluster.
    \f
    :param click.core.Context ctx: Click context.
    """

    # get context parameters
    ctx.obj.update(ctx.params)
    ctx.obj["mode"] = bs_enums.input_parameters.RUN_MODE.DEREPLICATE

    # set start time and label
    set_start(ctx.obj)

    # workflow validations
    validate_output_paths(ctx)

    # initiate logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # parse config file
    logging.info("Using config file %s", ctx.obj["config_file_path"])
    BigscapeConfig.parse_config(ctx.obj["config_file_path"], ctx.obj["log_path"])

    # run BiG-SCAPE dereplicate
    run_bigscape_dereplicate(ctx.obj)
