""" Click parameters for the BiG-SCAPE Cluster CLI command """

# from python
import click

# from other modules
from big_scape.main import run_bigscape
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_binning_cluster_workflow,
    validate_pfam_path,
    set_start,
)


@click.command()
@common_all
@common_cluster_query
# comparison parameters
@click.option(
    "--hybrids_off",
    is_flag=True,
    help=(
        "Toggle to add BGCs with hybrid predicted classes/categories to each "
        "subclass instead of a hybrid class/network (e.g. a 'terpene-nrps' BGC "
        "would be added to both the terpene and NRPS classes/networks instead of "
        "the terpene.nrps network). "
        "Only works if --classify/--legacy_classify is selected."
    ),
)
# binning parameters
@click.option("--no_mix", is_flag=True, help=("Don't run the all-vs-all analysis."))
# networking parameters
@click.option(
    "--include_singletons",
    is_flag=True,
    help=("Include singletons in the network. (default: False)"),
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
    ctx.obj["query_bgc_path"] = None
    ctx.obj["mode"] = "Cluster"

    # workflow validations
    validate_binning_cluster_workflow(ctx)
    validate_pfam_path(ctx)
    validate_output_paths(ctx)

    # set start time and run label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE cluster
    run_bigscape(ctx.obj)
