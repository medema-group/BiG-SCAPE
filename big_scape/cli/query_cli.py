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
    validate_classify,
    validate_output_paths,
    validate_disk_only,
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
    "-q",
    "--query-bgc-path",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    required=True,
    callback=validate_query_bgc,
    help=(
        "Path to query BGC .gbk file. BiG-SCAPE will compare "
        "all BGCs records in the input and reference folders to the query"
        " in a one-vs-all mode."
    ),
)
@click.option(
    "-n",
    "--query-record-number",
    type=int,
    required=False,
    help=(
        "Query BGC record number. Used to select the specific record "
        "from the query BGC .gbk, and is only relevant when running "
        "--record-type cand_cluster, protocluster or proto_core."
        " Warning: if interleaved or chemical hybrid proto cluster/cores "
        "are merged (see config.yml), the relevant number is that of the "
        "first record of the merged cluster (the one with the lowest number). "
        "e.g. if records 1 and 2 get merged, the relevant number is 1. "
    ),
)
@click.option(
    "--propagate",
    is_flag=True,
    help=(
        "By default, BiG-SCAPE will only generate edges between the query and reference"
        " BGC records. With the propagate flag, BiG-SCAPE will go through multiple cycles of "
        "edge generation until no new reference BGCs are connected to the query "
        "connected component. For more details, see the Wiki."
    ),
)
@click.option(
    "--classify",
    type=click.Choice(["none", "class", "category"]),
    default="none",
    callback=validate_classify,
    help=(
        "By default BiG-SCAPE will compare the query BGC record against any other "
        "supplied reference BGC records regardless of antiSMASH  product class/category. "
        "Instead, select 'class' or 'category' to run analyses on one class-specific bin, "
        "in which case only reference BGC records with the same class/category as the "
        "query record will be compared. Can be used in combination with --legacy-weights "
        "for .gbks produced by antiSMASH version 6 or higher. For older antiSMASH versions "
        "or if --legacy-weights is not selected, BiG-SCAPE will use the generic 'mix' weights: "
        "{JC: 0.2, AI: 0.05, DSS: 0.75, Anchor boost: 2.0}. (default: none)"
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
    ctx.obj["mix"] = None
    ctx.obj["hybrids_off"] = False
    ctx.obj["mode"] = "Query"
    ctx.obj["exclude_categories"] = None
    ctx.obj["include_categories"] = None
    ctx.obj["exclude_classes"] = None
    ctx.obj["include_classes"] = None

    # workflow validations
    validate_pfam_path(ctx)
    validate_output_paths(ctx)
    validate_binning_query_workflow(ctx)
    validate_query_record(ctx)
    validate_disk_only(ctx)

    # set start time and label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE
    run_bigscape(ctx.obj)
