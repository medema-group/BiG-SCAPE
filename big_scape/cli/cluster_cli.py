""" Click parameters for the BiG-SCAPE Cluster CLI command """

# from python
import click

# from other modules
from big_scape.run_bigscape import run_bigscape
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all, common_cluster_query
from .cli_validations import (
    validate_output_paths,
    validate_disk_only,
    validate_binning_cluster_workflow,
    validate_pfam_path,
    validate_domain_include_list,
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
@click.option(
    "--mix",
    is_flag=True,
    help=(
        "Calculate distances using a 'mix' bin, wherein no classification is applied. "
        "This will do an all-vs-all comparison, and is likely going to take a long time. "
        "This bin will use weights from the 'mix' weights distribution. Warning: these "
        "weights are not recommended for use with the record types protocluster/protocore, "
        "as they have been optimized and validated only for the 'region' record type."
    ),
)
@click.option(
    "--legacy_classify",
    is_flag=True,
    help=(
        "Does not use antiSMASH BGC classes to run analyses on "
        "class-based bins, instead it uses BiG-SCAPE v1 predefined groups: "
        "PKS1, PKSOther, NRPS, NRPS-PKS-hybrid, RiPP, Saccharide, Terpene, Others. "
        "Will also use BiG-SCAPE v1 legacy_weights for distance calculations. "
        "This feature is available for backwards compatibility with "
        "antiSMASH versions up to v7. For higher antiSMASH versions, use "
        "at your own risk, as BGC classes may have changed. All antiSMASH "
        "classes that this legacy mode does not recognize will be grouped in "
        "'others'."
    ),
)
# networking parameters
@click.option(
    "--include_singletons",
    is_flag=True,
    help=("Include singletons in the network."),
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
    ctx.obj["propagate"] = True  # compatibility with query wrt cc generation
    ctx.obj["mode"] = "Cluster"

    # workflow validations
    validate_binning_cluster_workflow(ctx)
    validate_pfam_path(ctx)
    validate_domain_include_list(ctx)
    validate_output_paths(ctx)
    validate_disk_only(ctx)

    # set start time and run label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE cluster
    run_bigscape(ctx.obj)
