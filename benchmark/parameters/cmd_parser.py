"""Module contains command parser for parsing benchmark workflow arguments"""

# from python
from argparse import ArgumentParser


def parse_cmd(args):
    """Parse arguments from command line

    Arguments:
        args (list): list of command line arguments

    Returns:
        ArgumentsParser: object storing input arguments
    """
    parser = ArgumentParser(
        prog="BiG-SCAPE benchmark",
        description="Benchmarking the Biosynthetic Gene Similarity Clustering and Prospecting Engine",
    )

    parser.add_argument(
        "-i",
        "--computed_gcfs",
        dest="computed_gcfs",
        default=None,
        type=str,
        help="Path pointing to the ouput folder of a BiG-SCAPE run",
    )

    parser.add_argument(
        "-g",
        "--curated_gcfs",
        dest="curated_gcfs",
        default=None,
        type=str,
        help="Path pointing to file with curated GCF assignments",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        type=str,
        help="Output directory for benchmarking results",
    )

    return parser.parse_args(args)
