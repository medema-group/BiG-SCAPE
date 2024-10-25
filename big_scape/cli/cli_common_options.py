"""Click parameter that are shared between multiple BiG-SCAPE modules"""

# from python
import click
from multiprocessing import cpu_count
from pathlib import Path

import big_scape.paths as bs_paths

# from this module
from .cli_validations import (
    validate_profiling,
    validate_not_empty_dir,
    validate_input_mode,
    validate_alignment_mode,
    validate_extend_strategy,
    validate_includelist_all,
    validate_includelist_any,
    validate_gcf_cutoffs,
    validate_filter_gbk,
    validate_record_type,
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
        click.option(
            "--config_file_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            default=Path(bs_paths.DEFAULT_CONFIG_FILE),
            help="Path to BiG-SCAPE config file, which stores values for a "
            "series of advanced use parameters. (default: bundled big_scape/config.yml).",
        ),
        # diagnostic parameters
        click.option(
            "-v",
            "--verbose",
            is_flag=True,
            help=(
                "Prints more detailed information of each step in the analysis, "
                "output all kinds of logs, including debugging log info, and writes to logfile. "
                "Toggle to activate."
            ),
        ),
        click.option(
            "--quiet",
            is_flag=True,
            help="Don't print any log info to output, only write to logfile.",
        ),
        # run parameters
        click.option(
            "-l",
            "--label",
            default=None,
            type=str,
            help="A run label to be added to the output results folder name, as well as "
            "dropdown menu in the visualization page. "
            "By default, BiG-SCAPE runs will have a name such as YYYY-MM-DD_HH-MM-SS_[label]",
        ),
        click.option(
            "-c",
            "--cores",
            default=cpu_count(),
            type=int,
            help=(
                "Set the max number of cores available (default: use all available cores)."
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
            help="Path to output log file. (default: output_dir/timestamp.log).",
        ),
        click.option(
            "--no-db-dump",
            type=bool,
            is_flag=True,
            default=False,
            help="Do not dump the sqlite database to disk. This will speed up your run, "
            "but in case of a crashed run no info will be stored and you'll have to "
            "re-start the run from scratch",
        ),
        click.option(
            "--disk-only",
            type=bool,
            is_flag=True,
            default=False,
            help="Do not store any results in memory, only on disk. This is almost certainly "
            "slower than the default behaviour, but can be useful for very large runs or "
            "runs with limited memory.",
        ),
        click.option(
            "--db-only-output",
            type=bool,
            is_flag=True,
            default=False,
            help="Do not generate any output besides the data stored in the database. "
            "Suitable for advanced users that wish to only make use of the results "
            "stored in the SQLite database.",
        ),
        click.option(
            "--no-trees",
            type=bool,
            is_flag=True,
            default=False,
            help="Do not generate any GCF newick trees. Suitable for users that do not "
            "utilize our output visualization, but only make use of the output stored "
            "in the tsv files and/or SQLite database.",
        ),
        click.option(
            "--force-gbk",
            type=bool,
            is_flag=True,
            default=False,
            help="If GBK files are found without antiSMASH annotations, this adds a region covering "
            "the full sequence, and sets its product to 'other'. Warning: BiG-SCAPE still "
            "needs CDS features and a sequence feature to work with non-antiSMASH gbks. "
            "Furthermore, this feature is still under development, use at own risk.",
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
            help="Run profiler and output profile report. Note: currently only available for Linux systems.",
        ),
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
                "Where to look for input GBK files. "
                "recursive: search for gbk files recursively in input directory. "
                "flat: search for gbk files in input directory only. "
                "(default: recursive)."
            ),
        ),
        # TODO: adjust choices
        click.option(
            "-m",
            "--mibig_version",
            type=str,
            required=False,
            help="MIBiG release number (from 3.1 onwards). If not provided, MIBiG will not be "
            "included in the analysis. If required, BiG-SCAPE will download the "
            "MIBiG database to ./big_scape/MIBiG/mibig_antismash_<version>_gbk. "
            "(Advanced) Any custom MIBiG collection can be used as long as the expected "
            "folder is present.",
        ),
        click.option(
            "-r",
            "--reference_dir",
            callback=validate_not_empty_dir,
            type=click.Path(
                exists=True, file_okay=False, dir_okay=True, path_type=Path
            ),
            help="Path to directory containing user defined, non-MIBiG, antiSMASH processed reference BGCs.",
        ),
        click.option(
            "--include_gbk",
            type=str,
            default="cluster,region",
            callback=validate_filter_gbk,
            help=(
                "A comma separated list of strings. Only gbk files that have "
                "the string(s) in their filename will be used for the analysis "
                "(default: 'cluster,region'). Use an asterisk to accept every "
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
            "-p",
            "--pfam_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help="Path to Pfam database file.",
        ),
        click.option(
            # TODO: implement
            "--domain_includelist_all_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            callback=validate_includelist_all,
            help=(
                "Path to txt file with Pfam accessions. Only BGCs containing all "
                "the listed accessions will be analysed. In this file, each "
                "line contains a single Pfam accession (with an optional comment,"
                " separated by a tab). Lines starting with '#' are ignored. Pfam "
                "accessions are case-sensitive."
            ),
        ),
        click.option(
            # TODO: implement
            "--domain_includelist_any_path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            callback=validate_includelist_any,
            help=(
                "Path to txt file with Pfam accessions. Only BGCs containing any of "
                "the listed accessions will be analysed. In this file, each "
                "line contains a single Pfam accession (with an optional comment,"
                " separated by a tab). Lines starting with '#' are ignored. Pfam "
                "accessions are case-sensitive."
            ),
        ),
        click.option(
            "--legacy_weights",
            is_flag=True,
            help=(
                "Use BiG-SCAPE v1 class-based weights in distance calculations. "
                "If not selected, the distance metric will be based on the 'mix' "
                "weights distribution. Warning: these legacy weights are not recommended "
                "for use with the record types 'protocluster'/'protocore', as they have "
                "been optimized and validated only for the 'region' record type."
            ),
        ),
        click.option(
            "--alignment_mode",
            type=click.Choice(["global", "glocal", "local", "auto"]),
            default="glocal",
            callback=validate_alignment_mode,
            help=(
                "Alignment mode for each pair of gene clusters. 'global': the whole "
                "list of domains of each BGC are compared; 'local': Longest Common "
                "Subcluster mode. Redefine the subset of the domains used to "
                "calculate distance by trying to find the longest slice of common "
                "domain content per gene in both BGCs, then extend each slice. "
                "'glocal': Similar to local, but extension assumes full extension "
                "of the shortest upstream/downstream arms in a compared pair. "
                "'auto': use glocal when at least one of the BGCs in each pair "
                "has the 'contig_edge' annotation from antiSMASH v4+, otherwise "
                "use global mode on that pair. For an in depth description, see the wiki."
                " (default: glocal)."
            ),
        ),
        click.option(
            "--extend_strategy",
            type=click.Choice(["legacy", "greedy", "simple_match"]),
            default="legacy",
            callback=validate_extend_strategy,
            help="Strategy to extend BGCs. 'legacy' will use the original BiG-SCAPE extension strategy, "
            "while 'greedy' or 'simple_match' will use new extension strategies. For an in depth description,"
            " see the wiki. (default: legacy).",
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
                "0.25,0.5,1.0. (default: 0.3)."
            ),
        ),
        # output parameters
        click.option(
            "--profile_path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to output profile file. (default: output_dir/).",
        ),
        click.option(
            "-db",
            "--db_path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to sqlite db output file. (default: output_dir/output_dir.db).",
        ),
        # TODO: implement cand_cluster here and LCS-ext
        click.option(
            # TODO: double check that cand_cluster is proper implemented
            "--record_type",
            type=click.Choice(["region", "cand_cluster", "protocluster", "proto_core"]),
            default="region",
            callback=validate_record_type,
            help="Use a specific type of antiSMASH record for comparison. (default: region).",
        ),
    ]

    for opt in options[::-1]:
        fn = opt(fn)
    return fn
