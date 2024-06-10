""" Click parameters for the BiG-SCAPE Query CLI command """

# from python
import click
from pathlib import Path

# from other modules
from big_scape.run_bigscape import run_bigscape
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_query_bgc,
    validate_pfam_path,
    set_start,
    validate_binning_query_workflow,
    validate_query_record,
)


@click.command()
@common_all
@common_cluster_query
@click.option(
    "-q" "--query_bgc_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    required=True,
    callback=validate_query_bgc,
    help=(
        "Path to query BGC file. BiG-SCAPE will compare "
        "all BGCs in the input and reference folders to the query"
        " in a one-vs-all mode."
    ),
)
@click.option(
    "-n" "--query_record_number",
    type=int,
    required=False,
    help=(
        "Query BGC record number. Used to select the specific record "
        "from the query BGC gbk. Warning: interleaved or chemical hybrid proto "
        "cluster/cores are merged, and the relevant number is that of the "
        "first record of the merged cluster (the one with the lowest number). "
        "e.g. if records 1 and 2 get merged, the relevant number is 1. "
    ),
)
@click.option(
    "--skip_propagation",
    is_flag=True,
    help=(
        "Only generate edges between the query and reference BGCs. If not set, "
        "BiG-SCAPE will also propagate edge generation to reference BGCs. "
        "Warning: if the database already contains all edges, this will not work, "
        "and the output will still showcase all edges between nodes "
        "in the query connected component."
    ),
)
@click.pass_context
def query(ctx, *args, **kwarg):
    """
    BiG-SCAPE - QUERY

    Query BGC mode - BiG-SCAPE queries a set of BGCs
    based on a single BGC query in a one-vs-all comparison.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    # get context parameters
    ctx.obj.update(ctx.params)
    ctx.obj["no_mix"] = None
    ctx.obj["hybrids_off"] = False
    ctx.obj["mode"] = "Query"

    # workflow validations
    validate_pfam_path(ctx)
    validate_output_paths(ctx)
    validate_binning_query_workflow(ctx)
    validate_query_record(ctx)

    # set start time and label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE
    run_bigscape(ctx.obj)
