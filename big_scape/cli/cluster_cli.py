""" Click parameters for the BiG-SCAPE Cluster CLI command """

# from python
import click
from pathlib import Path

# from other modules
from big_scape.cluster import run_bigscape_cluster, run_bigscape
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_binning_workflow,
    validate_skip_hmmscan,
    validate_pfam_path,
    set_start,
)


@click.command()
@common_all
@common_cluster_query
# TODO: delete once query bgc mode is ready
@click.option(
    "--query_bgc_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    help=(
        "Path to query BGC file. BiG-SCAPE will compare, "
        "all input BGCs to the query in a one-vs-all mode."
    ),
)
# binning parameters
@click.option("--no_mix", is_flag=True, help=("Dont run the all-vs-all analysis"))
@click.option(
    "--legacy_classify",
    is_flag=True,
    help=(
        "Does not use antiSMASH/BGC classes to run analyses on "
        "class-based bins, instead it uses BiG-SCAPEv1 predefined groups."
    ),
)
@click.option(
    "--classify",
    is_flag=True,
    help=("Use antiSMASH/BGC classes to run analyses on " "class-based bins."),
)
# networking parameters
@click.option(
    "--include_singletons",
    is_flag=True,
    help=("Include singletons in the network. Default: False"),
)
@click.pass_context
def cluster(ctx, *args, **kwargs):
    """
    BiG-SCAPE - CLUSTER

    Clustering mode - BiG-SCAPE performs clustering of BGCs into GCFs.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    # get context parameters
    ctx.obj.update(ctx.params)

    # workflow validations
    validate_binning_workflow(ctx)
    validate_skip_hmmscan(ctx)
    validate_pfam_path(ctx)
    validate_output_paths(ctx)

    # set start time and run label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE cluster
    run_bigscape_cluster(ctx.obj)
    run_bigscape(ctx.obj)
