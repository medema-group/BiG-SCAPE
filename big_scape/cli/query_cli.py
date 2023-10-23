""" Click parameters for the BiG-SCAPE Query CLI command """

# from python
import click
from pathlib import Path

# from other modules
from big_scape.query import run_bigscape_query

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_skip_hmmscan,
    validate_query_bgc,
    validate_pfam_path,
    set_start,
)


@click.command()
@common_all
@common_cluster_query
@click.option(
    "--query_bgc_path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    required=True,
    callback=validate_query_bgc,
    help=(
        "Path to query BGC file. BiG-SCAPE will compare, "
        "all input BGCs to the query in a one-vs-all mode."
    ),
)
@click.pass_context
def query(ctx, *args, **kwarg):
    """
    BiG-SCAPE - QUERY

    Query BGC mode - BiG-SCAPE queries a set of BGCs
    based on a single BGC query.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    # get context parameters
    ctx.obj.update(ctx.params)

    # workflow validations
    validate_skip_hmmscan(ctx)
    validate_pfam_path(ctx)
    validate_output_paths(ctx)

    # set start time and label
    set_start(ctx.obj)

    # run BiG-SCAPE query
    run_bigscape_query(ctx.obj)
