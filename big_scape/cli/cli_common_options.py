"""Click parameter that are shared between multiple BiG-SCAPE modules"""

# from python
import click
from multiprocessing import cpu_count
from pathlib import Path

# from this module
from .cli_validations import (
    validate_profiling,
    validate_not_empty_dir,
    validate_input_mode,
    validate_alignment_mode,
    validate_includelist,
    validate_gcf_cutoffs,
    validate_filter_gbk,
    validate_record_type,
    validate_classify,
    validate_output_dir,
)


def common_all(fn):
    """Decorator incorporating common command line options used in all modules

    Args:
        fn (function): Function to be decorated

    Returns:
        function decorated with new click options
    """
    options = [
        # meta
        click.option(
            "--metadata_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help="Path to metadata file.",
        ),
        click.option(
            "--config_file_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help="Path to BiG-SCAPE config file.",
        ),
        # diagnostic parameters
        click.option(
            "-v",
            "--verbose",
            is_flag=True,
            help=(
                "output all kinds of logs, "
                "including debugging log info, and write to logfile."
            ),
        ),
        click.option(
            "--quiet",
            is_flag=True,
            help="Don't print any log info to output, only write to logfile.",
        ),
        # run parameters
        click.option(
            "--label",
            default=None,
            type=str,
            help="A run label to be added to the output results folder name.",
        ),
        click.option(
            "-c",
            "--cores",
            default=cpu_count(),
            type=int,
            help=(
                "Set the max number of cores available"
                " (default: use all available cores)."
            ),
        ),
        # output parameters
        click.option(
            "-o",
            "--output_dir",
            type=click.Path(path_type=Path, dir_okay=True, file_okay=False),
            required=True,
            callback=validate_output_dir,
            help="Output directory for all BiG-SCAPE results files.",
        ),
        click.option(
            "--log_path",
            type=click.Path(
                path_type=Path(exists=False), dir_okay=True, file_okay=False
            ),
            help="Path to output log file. Default: output_dir/timestamp.log.",
        ),
        click.option(
            "--no-db-dump",
            type=bool,
            is_flag=True,
            default=False,
            help="Do not dump the sqlite database to disk",
        ),
    ]
    for opt in options[::-1]:
        fn = opt(fn)
    return fn


def common_cluster_query(fn):
    """Decorator incorporating common command line options in cluster and query modules

    Args:
        fn (function): Function to be decorated

    Returns:
        function decorated with new click options
    """
    options = [
        # diagnostic parameters
        click.option(
            "--profiling",
            callback=validate_profiling,
            is_flag=True,
            help="Run profiler and output profile report",
        ),
        # input parameters
        click.option(
            "-i",
            "--input_dir",
            "--gbk_dir",
            callback=validate_not_empty_dir,
            type=click.Path(
                exists=True, file_okay=False, dir_okay=True, path_type=Path
            ),
            required=True,
            help="Input directory containing gbk files to be used by BiG-SCAPE.",
        ),
        click.option(
            "--input_mode",
            default="recursive",
            callback=validate_input_mode,
            type=click.Choice(["recursive", "flat"]),
            help=(
                "Where to look for input GBK files. Default: recursive"
                "recursive: search for gbk files recursively in input directory. "
                "flat: search for gbk files in input directory only."
            ),
        ),
        # TODO: adjust choices
        click.option(
            "--mibig_version",
            type=click.Choice(
                ["1.0", "1.1", "1.2", "1.3", "1.4", "2.0", "3.0", "3.1", "custom"]
            ),
            required=False,
            help="MIBiG release number (e.g. 3.1). Download (if needed) and use this "
            "version of MiBIG database. If not provided, MiBIG will not be included "
            "in the analysis. If download is required, BiG-SCAPE "
            "will download the MIBiG database to the BiG-SCAPE folder",
        ),
        click.option(
            "--reference_dir",
            callback=validate_not_empty_dir,
            type=click.Path(
                exists=True, file_okay=False, dir_okay=True, path_type=Path
            ),
            help="Path to directory containing antismash-processed reference BGCs.",
        ),
        # TODO: check if still needed
        click.option(
            "--dataset_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help="Path to location of input dataset mapping file.",
        ),
        click.option(
            "--include_gbk",
            type=str,
            default="cluster,region",
            callback=validate_filter_gbk,
            help=(
                "A comma separated list of strings. "
                "Only gbk files with this string(s) will be used for the analysis "
                "(default: 'cluster', 'region'). Use an asterisk to accept every "
                "file (overrides '--exclude_gbk_str')."
            ),
        ),
        click.option(
            "--exclude_gbk",
            type=str,
            default="final",
            callback=validate_filter_gbk,
            help=(
                "A comma separated list of strings. "
                "If any string in this list occurs in the gbk filename, this "
                "file will not be used for the analysis (default: final)."
            ),
        ),
        click.option(
            "--min_bgc_length",
            type=int,
            default=0,
            help="Minimum BGC length to be included in analysis (default: 0bp).",
        ),
        click.option(
            "--cds_overlap_cutoff",
            type=click.FloatRange(min=0, max=1),
            default=0.1,
            help=(
                "Specify at which overlap percentage (as a decimal) two CDS in a gbk"
                " are considered to overlap. This preserves longest overlapping CDS."
            ),
        ),
        click.option(
            "-p",
            "--pfam_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help="Path to Pfam database file(s).",
        ),
        # hmmer parameters
        click.option(
            "--domain_overlap_cutoff",
            type=click.FloatRange(min=0, max=1),
            default=0.1,
            help=(
                "Specify at which overlap percentage (as a decimal) two domains"
                " in a CDS are considered to overlap. Domain with the "
                "best score is kept (default=0.1)"
            ),
        ),
        click.option(
            "--force_hmmscan",
            is_flag=True,
            help=(
                "Force domain prediction using hmmscan even if BiG-SCAPE finds "
                "processed gbk files (e.g. to use a new version of Pfam)."
            ),
        ),
        click.option(
            "--skip_hmmscan",
            is_flag=True,
            help=(
                "Skip domain prediction using hmmscan."
                "BiG-SCAPE expects to find "
                "a database of already processed gbks."
            ),
        ),
        click.option(
            "--domain_includelist_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            callback=validate_includelist,
            help=(
                "Path to txt file with Pfam accessions. Only BGCs containing "
                "the listed accessions will be analysed."
            ),
        ),
        click.option(
            "--legacy_weights",
            is_flag=True,
            help=(
                "Use BiG-SCAPE v1 class-based weights in distance calculations"
                "If not selected, the distance metric will be based on the 'mix'"
                " weights distribution."
            ),
        ),
        click.option(
            "--classify",
            type=click.Choice(["class", "category"]),
            callback=validate_classify,
            help=(
                "Use antiSMASH/BGC classes or categories to run analyses on class-based bins."
                "Can be used in combination with --legacy_weights if BGC gbks "
                "have been produced by antiSMASH version6 or higher. For older "
                "antiSMASH versions, either use --legacy_classify or do not select"
                "--legacy_weights, which will perform the weighted distance calculations"
                "based on the generic 'mix' weights."
            ),
        ),
        click.option(
            "--hybrids_off",
            is_flag=True,
            help=(
                "Toggle to add BGCs with hybrid predicted classes/categories to each"
                "subclass instead of a hybrid class/network (e.g. a 'terpene-nrps' BGC "
                "would be added to both the terpene and NRPS classes/networks instead of"
                " the terpene.nrps network)."
                " Only works if --classify/--legacy_classify is selected."
            ),
        ),
        click.option(
            "--alignment_mode",
            type=click.Choice(["global", "glocal", "auto"]),
            # TODO: define proper default
            default="auto",
            callback=validate_alignment_mode,
            help=(
                "Alignment mode for each pair of gene clusters. 'global': the whole "
                "list of domains of each BGC are compared; 'glocal': Longest Common "
                "Subcluster mode. Redefine the subset of the domains used to "
                "calculate distance by trying to find the longest slice of common "
                "domain content per gene in both BGCs, then expand each slice."
                " 'auto': use glocal when at least one of the BGCs in each pair "
                "has the 'contig_edge'annotation from antiSMASH v4+, otherwise "
                "use global mode on that pair"
            ),
        ),
        # networking parameters
        click.option(
            "--gcf_cutoffs",
            type=str,
            default=0.3,
            callback=validate_gcf_cutoffs,
            help=(
                "A comma separated list of floats. "
                "Generate networks using multiple raw distance cutoff values. "
                "Values should be in the range [0.0, 1.0]. Example: --gcf_cutoffs 0.1,"
                "0.25,0.5,1.0. Default: 0.3"
            ),
        ),
        # output parameters
        click.option(
            "--profile_path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to output profile file. Default: output_dir/.",
        ),
        click.option(
            "--db_path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to sqlite db output file. Default: output_dir/data_sqlite.db.",
        ),
        click.option(
            "--record_type",
            type=click.Choice(["region", "proto_cluster", "proto_core"]),
            default="region",
            callback=validate_record_type,
            help="Use a specific type of record for comparison. Default: region",
        ),
    ]
    for opt in options[::-1]:
        fn = opt(fn)
    return fn
