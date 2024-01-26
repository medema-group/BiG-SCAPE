#!/usr/bin/env python3

"""
BiG-SCAPE

PI: Marnix Medema                   marnix.medema@wur.nl

Maintainers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Loureiro                   catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl

Developers:
Jorge Navarro                       jorge.navarromunoz@wur.nl
Satria Kautsar                      sakautsar@lbl.gov
Emmanuel (Emzo) de los Santos       E.De-Los-Santos@warwick.ac.uk
Arjan Draisma                       arjan.draisma@wur.nl
Catarina Loureiro                   catarina.salesesantosloureiro@wur.nl
Nico Louwen                         nico.louwen@wur.nl


Usage:   Please see `python bigscape.py -h`
                    `python bigscape.py cluster -h`

Example: python bigscape.py cluster -i ./inputfiles -o ./results -c 8 -p ./Pfam-A.hmm


Official repository:
https://github.com/medema-group/BiG-SCAPE


# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in
# LICENSE.txt.
"""

# from python
import click

# from other modules
from big_scape.cli import cluster_cli as cluster
from big_scape.cli import query_cli as query
from big_scape.cli import benchmark_cli as benchmark


def print_version(ctx, param, value):
    if not value or ctx.resilient_parsing:
        return
    click.echo("BiG-SCAPE 2.0 beta")
    ctx.exit()


@click.group(context_settings={"help_option_names": ["-h", "--help"]})
@click.option(
    "--version", is_flag=True, callback=print_version, expose_value=False, is_eager=True
)
@click.pass_context
def cli(ctx, *args, **kwargs):
    """
    BiG-SCAPE

    Biosynthetic Gene Similarity Clustering and Prospecting Engine.

    BiG-SCAPE can be run in three modes: cluster, query, and benchmark.
    See the help menu for each mode for more information.
    For a more comprehensive help menu and tutorials see GitHub Wiki.
    \f
    :param click.core.Context ctx: Click context.
    """
    ctx.obj = ctx.params
    pass


cli.add_command(cluster.cluster)
cli.add_command(query.query)
cli.add_command(benchmark.benchmark)


if __name__ == "__main__":
    cli()
