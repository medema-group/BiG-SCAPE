"""Module which contains the cmd parser used in parsing and storing command
line arguments"""

# from python
from argparse import ArgumentParser
from multiprocessing import cpu_count
from pathlib import Path


def cmd_parser():
    """
    Parse arguments from the command line

    Returns:
        ArgumentParser: object storing input arguments
    """

    parser = ArgumentParser(
        prog="BiG-SCAPE",
        description="Biosynthetic Gene Similarity Clustering and Prospecting Engine",
        epilog="For a more comprehensive help menu and tutorials see GitHub Wiki",
    )

    parser.add_argument(
        "--label",
        dest="label",
        default=None,
        type=str,
        help="a run label to be added to the output results folder name."
        # TODO: check this when refactoring method
    )

    parser.add_argument(
        "-c",
        "--cores",
        dest="cores",
        default=cpu_count(),
        type=int,
        help="Set the max number of cores available (default: use all available\
        cores).",
    )

    parser.add_argument(
        "-i",
        "--inputdir",
        dest="inputdir",
        type=Path,
        required=True,
        help="Input directory containing gbk files to be used by BiG-SCAPE.",
    )

    parser.add_argument(
        "--input_mode",
        dest="input_mode",
        default="recursive",
        type=str,
        choices=["flat", "recursive"],
        help="Where to look for input GBK files. Default: recursive",
    )

    parser.add_argument(
        "-o",
        "--outputdir",
        dest="outputdir",
        type=Path,
        required=True,
        help="Output directory for all BiG-SCAPE results files.",
    )

    parser.add_argument(
        "--db_path",
        dest="db_path",
        type=Path,
        help="Path to sqlite db output directory. Default: outputdir.",
    )

    parser.add_argument(
        "--log_path",
        dest="log_path",
        type=Path,
        help="Path to output log file directory. Default: outputdir.",
    )

    parser.add_argument(
        "--metadata_path",
        dest="metadata_path",
        type=Path,
        help="Path to location of input metadata file.",
    )

    parser.add_argument(
        "--dataset_path",
        dest="dataset_path",
        type=Path,
        help="Path to location of input dataset mapping file."
        # TODO: check this when implementing method
    )

    parser.add_argument(
        "-p",
        "--pfam_path",
        dest="pfam_path",
        type=Path,
        required=True,
        help="Path to hmmpress-processed Pfam database files."
        # TODO: check this when implementing automated pressing
    )

    parser.add_argument(
        "--download_pfam",
        dest="download_pfam",
        default=False,
        action="store_true",
        help="Download and press latest pfam database.",
    )

    parser.add_argument(
        "--mibig_version",
        dest="mibig_version",
        type=str,
        help="Version of MiBIG database. If not provided, MiBIG will not be included\
             in the analysis.",
    )

    parser.add_argument(
        "--reference_dir",
        dest="reference_dir",
        type=Path,
        help="Directory containing user provided reference BGC antiSMASH-processed \
            genbank files.",
    )

    parser.add_argument(
        "--include_gbk",
        dest="include_gbk",
        default=["cluster", "region"],
        type=str,
        nargs="+",
        help="Only gbk files with this string(s) will be used for the analysis \
            (default: 'cluster', 'region'). Use an asterisk to accept every file \
            (overrides '--exclude_gbk_str').",
    )

    parser.add_argument(
        "--exclude_gbk",
        dest="exclude_gbk",
        default=["final"],
        type=str,
        nargs="+",
        help="If any string in this list occurs in the gbk filename, this \
            file will not be used for the analysis (default: final).",
    )

    parser.add_argument(
        "--query_bgc_path",
        dest="query_bgc_path",
        type=Path,
        help="Path to location of query BGC gbk file. When provided, BiG-SCAPE\
            will compare, all input BGCs to the query in a one-vs-all mode.",
    )

    parser.add_argument(
        "--min_bgc_length",
        dest="min_bgc_length",
        default=0,
        type=int,
        help="Provide the minimum size of a BGC to be included in the analysis.\
              Default is 0 base pairs.",
    )

    parser.add_argument(
        "--domain_overlap_cutoff",
        dest="domain_overlap_cutoff",
        default=0.1,
        type=float,
        help="Specify at which overlap percentage domains are considered to \
            overlap. Domain with the best score is kept (default=0.1).",
    )

    parser.add_argument(
        "--force_hmmscan",
        dest="force_hmmscan",
        default=False,
        action="store_true",
        help="Force domain prediction using hmmscan even if BiG-SCAPE finds \
            processed domtable files (e.g. to use a new version of Pfam).",
    )

    parser.add_argument(
        "--force_hmmalign",
        dest="force_hmmalign",
        default=False,
        action="store_true",
        help="Force multiple alignment of domains' sequences, even if alignments \
            have been generated in a previous run.",
    )

    parser.add_argument(
        "--domain_includelist_path",
        dest="domain_includelist_path",
        type=Path,
        help="Path to txt file with Pfam accessions. Only BGCs containing the \
            listed accessions will be analysed.",
    )

    parser.add_argument(
        "--mix",
        dest="mix",
        default=False,
        action="store_true",
        help="Run an all-vs-all analysis"
        # TODO: update with binning modes
    )

    parser.add_argument(
        "--alignment_mode",
        dest="alignment_mode",
        default="glocal",
        choices=["global", "glocal", "auto"],
        help="Alignment mode for each pair of gene clusters. 'global': the whole \
            list of domains of each BGC are compared; 'glocal': Longest Common \
            Subcluster mode. Redefine the subset of the domains used to calculate \
            distance by trying to find the longest slice of common domain content \
            per gene in both BGCs, then expand each slice. 'auto': use glocal when \
            at least one of the BGCs in each pair has the 'contig_edge' annotation \
            from antiSMASH v4+, otherwise use global mode on that pair"
        # TODO: update with implementation
    )

    parser.add_argument(
        "--gcf_cutoffs",
        dest="gcf_cutoffs",
        default="0.30",
        type=str,
        help="Generate networks using multiple raw distance cutoff values. \
            Values should be in the range [0.0, 1.0]. Example: --cutoffs 0.1,\
            0.25,0.5,1.0. Default: 0.3.",
    )

    parser.add_argument(
        "--include_singletons",
        dest="include_singletons",
        default=False,
        action="store_true",
        help="Include nodes that have no edges to other nodes from the network.",
    )

    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        default=False,
        action="store_true",
        help="output all kinds of logs, including debugging log info, and write\
            to logfile.",
    )

    parser.add_argument(
        "--quiet",
        dest="quiet",
        default=False,
        action="store_true",
        help="Don't print any log info to output, only write to logfile.",
    )

    parser.add_argument(
        "--profiling",
        dest="profiling",
        default=False,
        action="store_true",
        help="Run profiler and output profile report",
    )

    parser.add_argument(
        "--version",
        dest="version",
        action="version",
        version="%(prog)s 2.0 beta",
        help="Output software version info "
        # TODO: update with release
    )

    return parser
