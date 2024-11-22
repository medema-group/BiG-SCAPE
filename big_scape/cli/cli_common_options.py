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
            "--config-file-path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            default=Path(bs_paths.DEFAULT_CONFIG_FILE),
            help=(
                "Path to BiG-SCAPE config.yml file, which stores values for a "
                "series of advanced use parameters. (default: bundled big_scape/config.yml)."
            ),
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
            help=(
                "A run label to be added to the output results folder name, as well as "
                "dropdown menu in the visualization page. "
                "By default, BiG-SCAPE runs will have a name such as [label]_YYYY-MM-DD_HH-MM-SS."
            ),
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
            "--output-dir",
            type=click.Path(path_type=Path, dir_okay=True, file_okay=False),
            required=True,
            callback=validate_output_dir,
            help="Output directory for all BiG-SCAPE results files.",
        ),
        click.option(
            "--log-path",
            type=click.Path(
                path_type=Path(exists=False), dir_okay=True, file_okay=False
            ),
            help="Path to output log file. (default: output_dir/timestamp.log).",
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
            "--input-dir",
            "--gbk-dir",
            callback=validate_not_empty_dir,
            type=click.Path(
                exists=True, file_okay=False, dir_okay=True, path_type=Path
            ),
            required=True,
            help="Input directory containing .gbk files to be used by BiG-SCAPE. See the wiki for more details.",
        ),
        click.option(
            "--input-mode",
            default="recursive",
            callback=validate_input_mode,
            type=click.Choice(["recursive", "flat"]),
            help=(
                "Tells BiG-SCAPE Where to look for input GBK files. "
                "recursive: search for .gbk files recursively in input directory. "
                "flat: search for .gbk files in input directory only. "
                "(default: recursive)."
            ),
        ),
        # TODO: adjust choices
        click.option(
            "-m",
            "--mibig-version",
            type=str,
            required=False,
            help=(
                "MIBiG release number (from 3.1 onwards). If not provided, MIBiG will not be "
                "included in the analysis. If required, BiG-SCAPE will download the "
                "MIBiG database to ./big_scape/MIBiG/mibig_antismash_<version>_gbk. "
                "For advanced users: any custom (antiSMASH-processed) MIBiG collection"
                " can be used as long as the expected folder is present, e.g. provide"
                " -m mymibig with ./big_scape/MIBiG/mibig_antismash_mymibig_gbk."
                "For more information, see the wiki."
            ),
        ),
        click.option(
            "-r",
            "--reference-dir",
            callback=validate_not_empty_dir,
            type=click.Path(
                exists=True, file_okay=False, dir_okay=True, path_type=Path
            ),
            help=(
                "Path to directory containing user defined, non-MIBiG, antiSMASH processed reference BGCs."
                "For more information, see the wiki."
            ),
        ),
        click.option(
            "--include-gbk",
            type=str,
            default="cluster,region",
            callback=validate_filter_gbk,
            help=(
                "A comma separated list of strings. Only .gbk files that have "
                "the string(s) in their filename will be used for the analysis."
                "Use an asterisk to accept every file ('*' overrides '--exclude_gbk_str')."
                "(default: cluster, region)."
            ),
        ),
        click.option(
            "--exclude-gbk",
            type=str,
            default="final",
            callback=validate_filter_gbk,
            help=(
                "A comma separated list of strings. "
                "If any string in this list occurs in the .gbk filename, this "
                "file will not be used for the analysis (default: final)."
            ),
        ),
        click.option(
            "-p",
            "--pfam-path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            help=(
                "Path to Pfam database `.hmm` file (e.g `Pfam-A.hmm`)."
                " If the `.hmm` file has already been pressed and the pressed files"
                " are included in the same folder as the Pfam `.hmm` file, BiG-SCAPE "
                "will also use these pressed files. If this is not the case, BiG-SCAPE"
                " will run `hmmpress`. Note: the latter requires the user to have write "
                "permissions to the given Pfam folder."
            ),
        ),
        click.option(
            # TODO: implement
            "--domain-includelist-all-path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            callback=validate_includelist_all,
            help=(
                "Path to .txt file with phmm domain accessions (commonly, Pfam accessions "
                "(e.g. PF00501)). Only regions containing all the listed accessions will "
                "be analyzed. In this file, each line contains a single phmm domain accession "
                "(with an optional comment, separated by a tab). Lines starting with '#' "
                "are ignored. Domain accessions are case-sensitive. Cannot be provided in "
                "conjuction with --domain-includelist-any-path."
            ),
        ),
        click.option(
            # TODO: implement
            "--domain-includelist-any-path",
            type=click.Path(
                exists=True, dir_okay=False, file_okay=True, path_type=Path
            ),
            callback=validate_includelist_any,
            help=(
                "Path to .txt file with phmm domain accessions (commonly, Pfam accessions "
                "(e.g. PF00501)). Only BGCs containing any of the listed accessions will "
                "be analyzed. In this file, each line contains a single phmm domain accession "
                "(with an optional comment, separated by a tab). Lines starting with '#' "
                "are ignored.  Domain accessions are case-sensitive. Cannot be provided in "
                "conjuction with --domain-includelist-all-path."
            ),
        ),
        click.option(
            "--legacy-weights",
            is_flag=True,
            help=(
                "Use BiG-SCAPE v1 class-based weights in distance calculations. "
                "If not selected, the distance metric will be based on the 'mix' "
                "weights distribution. Warning: these weights have not been validated "
                "for record types other than region (see option --record_type, and wiki)."
            ),
        ),
        click.option(
            "--alignment-mode",
            type=click.Choice(["global", "glocal", "local", "auto"]),
            default="glocal",
            callback=validate_alignment_mode,
            help=(
                "Alignment mode for each pair of gene clusters. global: the whole list of domains of "
                "each BGC record is compared; local: Seeds the subset of the domains used to calculate "
                "distance by trying to find the longest slice of common domain content (Longest Common "
                "Subcluster, LCS) between both records, then extends each side (see --extension_strategy). "
                "glocal: Starts with performing local, but domain selection is then extended to the "
                "shortest upstream/downstream arms in a compared record pair. "
                "auto: use glocal when at least one of the BGCs in each pair has the contig_edge annotation "
                "from antiSMASH v4+, otherwise use global mode on that pair. For an in depth description, see the wiki."
                " (default: glocal)."
            ),
        ),
        click.option(
            "--extend-strategy",
            type=click.Choice(["legacy", "greedy", "simple_match"]),
            default="legacy",
            callback=validate_extend_strategy,
            help=(
                "Strategy to extend the BGCs record pair comparable region. legacy will use the original "
                "BiG-SCAPE extend strategy, while greedy and simple match are newly introduced in BiG-SCAPE 2. "
                "Legacy and simple match both examine the domain architecture of the record pair in order to "
                "find the most relevant extended borders. Greedy is a very simple method that takes the "
                "coordinates of the outermost matching domains as the extended borders. "
                "For more detail see the wiki. (default: legacy)."
            ),
        ),
        # networking parameters
        click.option(
            "--gcf-cutoffs",
            type=str,
            default=0.3,
            callback=validate_gcf_cutoffs,
            help=(
                "A comma separated list of floats. "
                "Generate networks using multiple distance cutoff values. "
                "Values should be in the range [0.0, 1.0]. Example: --gcf_cutoffs 0.1,"
                "0.25,0.5,1.0. For more detail see the wiki. (default: 0.3)."
            ),
        ),
        # output parameters
        click.option(
            "--profile-path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to output profile file. (default: output_dir/).",
        ),
        click.option(
            "-db",
            "--db-path",
            type=click.Path(path_type=Path, dir_okay=False),
            help="Path to sqlite db output file. (default: output_dir/output_dir.db).",
        ),
        # TODO: implement cand_cluster here and LCS-ext
        click.option(
            # TODO: double check that cand_cluster is proper implemented
            "--record-type",
            type=click.Choice(["region", "cand_cluster", "protocluster", "proto_core"]),
            default="region",
            callback=validate_record_type,
            help=(
                "Use a specific type of antiSMASH record for comparison. For every .gbk, "
                "BiG-SCAPE will try to extract the requested record type, if this is not present, "
                "BiG-SCAPE will try to extract the next higher level record type, i.e. if a "
                "proto_core feature is not present, BiG-SCAPE will look for a protocluster feature, "
                "and so on and so forth. The record type hierarchy is: region>cand_cluster>protocluster>proto_core."
                ". For more detail, see the wiki. (default: region)."
            ),
        ),
        click.option(
            "--no-db-dump",
            type=bool,
            is_flag=True,
            default=False,
            help=(
                "Do not dump the sqlite database to disk. This will speed up your run, "
                "but in case of a crashed run no info will be stored and you'll have to "
                "re-start the run from scratch"
            ),
        ),
        click.option(
            "--disk-only",
            type=bool,
            is_flag=True,
            default=False,
            help=(
                "Do not store any results in memory, only on disk. This is almost certainly "
                "slower than the default behavior, but will save some memory and can therefore be "
                "useful for very large runs or runs with limited memory"
            ),
        ),
        click.option(
            "--db-only-output",
            type=bool,
            is_flag=True,
            default=False,
            help=(
                "Do not generate any output besides the data stored in the database. "
                "Suitable for advanced users that wish to only make use of the results "
                "stored in the SQLite database."
            ),
        ),
        click.option(
            "--no-trees",
            type=bool,
            is_flag=True,
            default=False,
            help=(
                "Do not generate any GCF newick trees. Suitable for users that do not "
                "utilize our output visualization, but only make use of the output stored "
                "in the .tsv files (which include the network files) and/or SQLite database."
            ),
        ),
        click.option(
            "--force-gbk",
            type=bool,
            is_flag=True,
            default=False,
            help=(
                "Recommended for advanced users only. Allows BiG-SCAPE to use non-antiSMASH "
                "processed .gbk files. If GBK files are found without antiSMASH annotations "
                "(specifically, BiG-SCAPE checks for the absence of a antiSMASH version feature), "
                "BiG-SCAPE will still read and parse these files, and will create internal gbk "
                "record objects, each of which will have a region feature covering the full sequence "
                "length and a product feature `other`. Warning: BiG-SCAPE still needs CDS features and "
                "a sequence feature to work with non-antiSMASH .gbks. Furthermore, --include-gbk and "
                "--exclude-gbk parameters might need to be adjusted if .gbk file names also do not follow "
                "antiSMASH format. Disclaimer: this feature is still under development, use at own risk."
            ),
        ),
    ]

    for opt in options[::-1]:
        fn = opt(fn)
    return fn
