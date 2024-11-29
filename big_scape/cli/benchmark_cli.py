""" Click parameters for the BiG-SCAPE Benchmark CLI command """

# from python
import click
from pathlib import Path

# from other modules
from big_scape.benchmark import run_bigscape_benchmark
from big_scape.diagnostics import init_logger, init_logger_file

# from this module
from .cli_common_options import common_all
from .cli_validations import (
    set_start,
    validate_output_paths,
)


# BiG-SCAPE benchmark mode
@click.command()
@common_all
# input parameters
@click.option(
    "-g",
    "--GCF-assignment-file",
    type=click.Path(exists=True, dir_okay=False, file_okay=True, path_type=Path),
    required=True,
    help=(
        "Path to GCF assignments file. BiG-SCAPE will compare "
        "a run output to these assignments."
    ),
)
@click.option(
    "-b",
    "--BiG-dir",
    type=click.Path(exists=True, dir_okay=True, file_okay=False, path_type=Path),
    required=True,
    help="Path to BiG-SCAPE (v1 or v2) or BiG-SLICE output directory.",
)
@click.pass_context
def benchmark(ctx, *args, **kwargs):
    """
    BiG-SCAPE - BENCHMARK

    Benchmarking mode - Compare a BiG-SCAPE or BiG-SLICE BGC clustering against
    a known/expected set of BGC <-> GCF assignments.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    # get context parameters
    ctx.obj.update(ctx.params)
    ctx.obj["mode"] = "Benchmark"

    # workflow validations
    validate_output_paths(ctx)

    # set start time and label
    set_start(ctx.obj)

    # initialize logger
    init_logger(ctx.obj)
    init_logger_file(ctx.obj)

    # run BiG-SCAPE benchmark
    run_bigscape_benchmark(ctx.obj)
