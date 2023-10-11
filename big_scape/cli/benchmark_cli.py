""" Click parameters for the BiG-SCAPE Benchmark CLI command """

# from python
import click
from multiprocessing import cpu_count
from pathlib import Path

# from other modules
from big_scape.benchmark import run_bigscape_benchmark


# BiG-SCAPE benchmark mode
@click.command()
# meta
@click.option(
    "--metadata_path",
    type=click.Path(exists=True, file_okay=True, path_type=Path),
    help="Path to metadata file.",
)
@click.option(
    "--config_file_path",
    type=click.Path(exists=True, file_okay=True, path_type=Path),
    help="Path to BiG-SCAPE config file.",
)
# diagnostic parameters
@click.option(
    "-v",
    "--verbose",
    is_flag=True,
    help=(
        "output all kinds of logs, "
        "including debugging log info, and write to logfile."
    ),
)
@click.option(
    "--quiet",
    is_flag=True,
    help="Don't print any log info to output, only write to logfile.",
)
# run parameters
@click.option(
    "--label",
    default=None,
    type=str,
    help="A run label to be added to the output results folder name.",
)
@click.option(
    "-c",
    "--cores",
    default=cpu_count(),
    type=int,
    help=(
        "Set the max number of cores available" " (default: use all available cores)."
    ),
)
# input parameters
@click.option(
    "--GCF_assigment_file",
    type=click.Path(exists=True, file_okay=True, path_type=Path),
    required=True,
    help=(
        "Path to GCF assignments file. BiG-SCAPE will compare, "
        "a run output to these assignments."
    ),
)
@click.option(
    "--run_database",
    type=click.Path(exists=True, file_okay=True, path_type=Path),
    required=True,
    help="Path to BiG-SCAPE run database.",
)
# output parameters
@click.option(
    "-o",
    "--output_dir",
    type=click.Path(path_type=Path),
    required=False,
    help="Output directory for all BiG-SCAPE results files.",
)
@click.option(
    "--log_path",
    type=click.Path(path_type=Path),
    help="Path to output log file directory. Default: output_dir.",
)
@click.pass_context
def benchmark(ctx, query, benchmark):
    """
    BiG-SCAPE - BENCHMARK

    Benchmarking mode - BiG-SCAPE compares the results of a query against
    a benchmark set of BGC <-> GCF assignments.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    ctx.obj.update(ctx.params)
    run_bigscape_benchmark(ctx.obj)
