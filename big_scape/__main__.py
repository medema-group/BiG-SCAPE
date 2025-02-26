"""Main file executed when the package is run as a module
"""

# from python
import click

# from other modules
from big_scape.cli import cluster_cli as cluster
from big_scape.cli import query_cli as query
from big_scape.cli import dereplicate_cli as dereplicate
from big_scape.cli import benchmark_cli as benchmark
from big_scape.cli.cli_config import CommandOrder
from big_scape.utility.version import get_bigscape_version


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    bigscape_version = get_bigscape_version()
    click.echo(f"BiG-SCAPE {bigscape_version}")
    ctx.exit()


@click.group(cls=CommandOrder, context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--version", is_flag=True, callback=print_version, expose_value=False, is_eager=True
)
@click.pass_context
def cli(ctx, *args, **kwargs):
    """
    BiG-SCAPE

    Biosynthetic Gene Similarity Clustering and Prospecting Engine.

    \b
    BiG-SCAPE can be run in four modes: cluster, query, dereplicate and benchmark.
    See the help menu for each mode for more information.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    ctx.obj = ctx.params
    pass


cli.add_command(cluster.cluster)
cli.add_command(query.query)
cli.add_command(dereplicate.dereplicate)
cli.add_command(benchmark.benchmark)


def main():
    cli()


if __name__ == "__main__":
    main()
