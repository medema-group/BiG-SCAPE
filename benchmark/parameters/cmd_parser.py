"""Contains command parser for parsing benchmark workflow arguments"""

# from python
from argparse import ArgumentParser
from pathlib import Path

# from other modules
from big_scape.errors import InvalidArgumentError


def parse_cmd(args):
    """Parse arguments from command line

    Args:
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
        "--input_dir",
        dest="input_dir",
        default=None,
        type=Path,
        required=True,
        help="Path pointing to the output folder of a BiG-SCAPE or BiG-SLICE run",
    )

    parser.add_argument(
        "-g",
        "--curated_gcfs",
        dest="curated_gcfs",
        default=None,
        type=Path,
        required=True,
        help="Path pointing to file with curated GCF assignments",
    )

    parser.add_argument(
        "-o",
        "--output_dir",
        dest="output_dir",
        type=Path,
        required=True,
        help="Output directory for benchmarking results",
    )

    return parser.parse_args(args)


def validate_args(args):
    """validate benchmarking arguments

    Args:
        ArgumentParser: object containing all required arguments

    Raises:
        InvalidArgumentError: upon encountering an invalid argument
    """
    if not args.input_dir.exists():
        raise InvalidArgumentError("--input_dir", args.input_dir)

    if not args.curated_gcfs.exists() or not args.curated_gcfs.is_file():
        raise InvalidArgumentError("--curated_gcfs", args.curated_gcfs)

    # empty file
    if args.curated_gcfs.stat().st_size == 0:
        raise InvalidArgumentError("--curated_gcfs", args.curated_gcfs)
