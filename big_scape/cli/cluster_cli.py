""" Click parameters for the BiG-SCAPE Cluster CLI command """

# from python
import click
from pathlib import Path

# from other modules
from big_scape.main import run_bigscape
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_binning_cluster_workflow,
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
# comparison parameters
@click.option(
    "--legacy_classify",
    is_flag=True,
    help=(
        "Does not use antiSMASH BGC classes to run analyses on "
        "class-based bins, instead it uses BiG-SCAPE v1 predefined groups: "
        "PKS1, PKSOther, NRPS, NRPS-PKS-hybrid, RiPP, Saccharide, Terpene, Others."
        "Will also use BiG-SCAPEv1 legacy_weights for distance calculations."
        "This feature is available for backwards compatibility with "
        "antiSMASH versions up to v7. For higher antiSMASH versions, use"
        " at your own risk, as BGC classes may have changed. All antiSMASH"
        "classes that this legacy mode does not recognize will be grouped in"
        " 'others'."
    ),
)
# binning parameters
@click.option("--no_mix", is_flag=True, help=("Dont run the all-vs-all analysis"))
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
    ctx.obj["query_bgc_path"] = None
    ctx.obj["mode"] = "Cluster"

    # workflow validations
    validate_binning_cluster_workflow(ctx)
    validate_skip_hmmscan(ctx)
    validate_pfam_path(ctx)
    validate_output_paths(ctx)

    # set start time and run label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE cluster
    print(ctx.obj)
    run_bigscape(ctx.obj)
