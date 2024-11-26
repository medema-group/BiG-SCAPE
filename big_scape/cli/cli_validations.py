""" A module to store all CLI parameter validations """

# from python
from datetime import datetime
import click
import logging
from pathlib import Path
from typing import Optional
import os
import time
import platform

# from other modules
from big_scape.errors.input_args import InvalidArgumentError

import big_scape.enums as bs_enums


# set time and label
# TODO: not a validation, move to another module
def set_start(param_dict) -> None:
    """get start time and set label in a run parameter dict"""

    start_time: datetime = datetime.now()
    timestamp = start_time.strftime("%d-%m-%Y_%H-%M-%S")
    if param_dict["label"]:
        param_dict["label"] = f"{param_dict['label']}_{timestamp}"
    else:
        param_dict["label"] = f"{timestamp}"

    param_dict["start_time"] = start_time

    return None


# meta parameter validations


def validate_profiling(ctx, param, profiling) -> bool:
    """Checks whether multithreading is possible, and therefore whether profiling can happen"""

    if profiling and platform.system() == "Darwin":
        logging.warning("Profiling is not supported on MacOS, please use Linux")
        return False
    return profiling


# input parameter validations


def validate_not_empty_dir(ctx, param, dir) -> Path:
    """Validates that a given directory is not empty.
    Raises a BadParameter"""

    if dir and dir.exists():
        contents = os.listdir(dir)
        if len(contents) == 0:
            logging.error(f"{dir}/ directory is empty!")
            raise click.BadParameter(f"{dir}/ directory is empty!")
    return dir


def validate_input_mode(ctx, param, input_mode) -> Optional[bs_enums.INPUT_MODE]:
    """validates the input_mode property. Raises an InvalidArgumentError if the
    input_mode parameter is invalid
    """
    # TODO: its not raising an error if the input_mode is invalid

    # check if the property matches one of the enum values
    valid_modes = [mode.value for mode in bs_enums.INPUT_MODE]

    for mode in valid_modes:
        if input_mode == mode:
            return bs_enums.INPUT_MODE[mode.upper()]
    return None


def validate_query_bgc(ctx, param, query_bgc_path) -> Path:
    """Raises an InvalidArgumentError if the query bgc path does not exist"""

    if query_bgc_path.suffix != ".gbk":
        logging.error(f"Query BGC file {query_bgc_path} is not a .gbk file!")
        raise click.BadParameter(f"Query BGC file {query_bgc_path} is not a .gbk file!")

    return query_bgc_path


def validate_class_category_filter(crx, param, filters) -> set[str]:
    """Parses provided include/exclude class/category filters into sets"""
    if filters is None:
        return set()
    return set(filters.split(","))


# output parameter validations


def validate_output_dir(ctx, param, output_dir) -> Path:
    """Validates that output directory exists"""

    if not output_dir.exists():
        parent_dir = output_dir.parent.absolute()
        if parent_dir.exists():
            os.makedirs(output_dir)
        else:
            logging.error(
                f"Output directory {output_dir} does not exist, and parent directory neither. Please create either."
            )
            raise click.BadParameter(
                f"Output directory {output_dir} does not exist, and parent directory neither. Please create either."
            )

    return output_dir


def validate_output_paths(ctx) -> None:
    """Sets the output paths to default output_dir if not provided"""

    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

    if "db_path" in ctx.obj and ctx.obj["db_path"] is None:
        db_path = ctx.obj["output_dir"] / Path(f"{ctx.obj['output_dir'].name}.db")
        ctx.obj["db_path"] = db_path

    if "log_path" in ctx.obj and ctx.obj["log_path"] is None:
        log_filename = timestamp + ".log"
        log_path = ctx.obj["output_dir"] / Path(log_filename)
        ctx.obj["log_path"] = log_path

    if "profile_path" in ctx.obj and ctx.obj["profile_path"] is None:
        profile_filename = timestamp + ".profile"
        profile_path = ctx.obj["output_dir"] / Path(profile_filename)
        ctx.obj["profile_path"] = profile_path

    return None


# db modes validations
def validate_disk_only(ctx) -> None:
    """Checks if the database storage/dumping modes that were set are compatible"""

    if not ("no_db_dump" in ctx.obj and "disk_only" in ctx.obj):
        raise RuntimeError(
            "Something went wrong with the database storage/dumping mode parameters. "
            "Please contact the developers."
        )

    if ctx.obj["no_db_dump"] and ctx.obj["disk_only"]:
        logging.error(
            "You have selected both --no-db-dump and --disk-only. Please select only one"
        )
        raise click.UsageError(
            "You have selected both --no-db-dump and --disk-only. Please select only one"
        )


# comparison validations


def validate_classify(ctx, param, classify) -> Optional[bs_enums.CLASSIFY_MODE]:
    """Validates whether the classification type is set, and if not
    sets the parameter to False"""

    # check if the property matches one of the enum values
    valid_modes = [mode.value for mode in bs_enums.CLASSIFY_MODE]

    if classify == "none":
        classify = False

    for mode in valid_modes:
        if classify == mode:
            return bs_enums.CLASSIFY_MODE[mode.upper()]

    return classify


def validate_alignment_mode(
    ctx, param, alignment_mode
) -> Optional[bs_enums.ALIGNMENT_MODE]:
    """Validate the passed alignment mode is one of the allowed modes"""
    # this should not ever run so long as the choices in the cmd_parser
    # remain updated and relevant
    valid_modes = [mode.value for mode in bs_enums.ALIGNMENT_MODE]

    for mode in valid_modes:
        if alignment_mode == mode:
            return bs_enums.ALIGNMENT_MODE[mode.upper()]
    return None


def validate_extend_strategy(
    ctx, param, extend_strategy
) -> Optional[bs_enums.EXTEND_STRATEGY]:
    """Validate the passed extend strategy is one of the allowed modes"""
    valid_strats = [strat.value for strat in bs_enums.EXTEND_STRATEGY]

    for strat in valid_strats:
        if extend_strategy == strat:
            return bs_enums.EXTEND_STRATEGY[strat.upper()]
    return None


def validate_gcf_cutoffs(ctx, param, gcf_cutoffs) -> list[float]:
    """Validates range and formats into correct list[float] format"""

    if "," in gcf_cutoffs:
        gcf_cutoffs = gcf_cutoffs.split(",")
        try:
            gcf_cutoffs = [float(cutoff) for cutoff in gcf_cutoffs]
        except ValueError:
            raise click.BadParameter("Invalid GCF cutoff value(s) provided!")
    else:
        try:
            gcf_cutoffs = [float(gcf_cutoffs)]
        except ValueError:
            raise click.BadParameter("Invalid GCF cutoff value(s) provided!")

    for cutoff in gcf_cutoffs:
        if cutoff < 0.0 or cutoff > 1.0:
            raise click.BadParameter("Invalid GCF cutoff value(s) provided!")

    return sorted(gcf_cutoffs, reverse=True)


def validate_filter_gbk(ctx, param, filter_str) -> list[str]:
    """Validates and formats the filter string and returns a list of strings"""

    return filter_str.split(",")


# hmmer parameters


def validate_includelist(ctx, param, domain_includelist_path):
    """Validate the path to the domain include list and return a list of domain
    accession strings contained within this file

    Returns:
        list[str]: A list of domain accessions to include
    """

    # only validate if set
    if domain_includelist_path is None:
        return None

    if not domain_includelist_path.exists():
        logging.error("domain-includelist file does not exist!")
        raise InvalidArgumentError(
            "--domain-includelist-all/any-path", domain_includelist_path
        )

    with domain_includelist_path.open(encoding="utf-8") as domain_includelist_file:
        lines = domain_includelist_file.readlines()

        pfams = []

        for line in lines:
            line = line.strip()
            elemts = line.split("\t")
            pfams.append(elemts[0])

        # expect Pfam accessions, i.e. PF00001 or PF00001.10
        lines_valid = map(
            lambda string: string.startswith("PF") and len(string) in range(7, 11),
            pfams,
        )

        if not all(lines_valid):
            logging.error(
                "Invalid Pfam accession(s) found in file %s", domain_includelist_path
            )
            raise click.BadParameter(
                "Invalid Pfam accession(s) found in file %s", domain_includelist_path
            )

        return pfams


def validate_includelist_all(ctx, param, domain_includelist_all_path) -> None:
    """Validate the path to the domain include list and return a list of domain
    accession strings contained within this file

    Returns:
        list[str]: A list of domain accessions to include
    """

    pfams = validate_includelist(ctx, param, domain_includelist_all_path)

    ctx.params["domain_includelist_all"] = pfams

    return None


def validate_includelist_any(ctx, param, domain_includelist_any_path) -> None:
    """Validate the path to the domain include list and return a list of domain
    accession strings contained within this file

    Returns:
        list[str]: A list of domain accessions to include
    """

    pfams = validate_includelist(ctx, param, domain_includelist_any_path)

    ctx.params["domain_includelist_any"] = pfams

    return None


# workflow validations


def validate_binning_cluster_workflow(ctx) -> None:
    """Raise an error if the combination of parameters in this object means no work
    will be done, or if the combination of parameters is invalid
    """

    # --legacy_weights needs a classification method

    if ctx.obj["legacy_weights"] and not ctx.obj["classify"]:
        logging.error(
            "You have selected --legacy-weights but no classification method. "
            "Please select any --classify method"
        )
        raise click.UsageError(
            "You have selected --legacy-weights but no classification method. "
            "Please select any --classify method"
        )

    # running without --mix needs a classification method, otherwise no work will be done

    if ctx.obj["mix"] is False and ctx.obj["classify"] is False:
        logging.error(
            "The combination of arguments you have selected for binning means no work will "
            "be done. Please run BiG-SCAPE using --mix, or add any --classify method"
            " in order to enable comparisons"
        )
        raise click.UsageError(
            "The combination of arguments you have selected for binning means no work will "
            "be done. Please run BiG-SCAPE using --mix, or add any --classify method"
            " in order to enable comparisons"
        )

    # --classify legacy turns on --legacy_weights

    if ctx.obj["classify"] == bs_enums.CLASSIFY_MODE.LEGACY:
        ctx.obj["legacy_weights"] = True

    # --hybrids_off needs --classify

    if ctx.obj["hybrids_off"]:
        if not (ctx.obj["classify"]):
            logging.error(
                "You have selected --hybrids-off but no classification method. "
                "Please select any --classify method"
            )
            raise click.UsageError(
                "You have selected --hybrids-off but no classification method. "
                "Please select any --classify method"
            )


def validate_binning_query_workflow(ctx) -> None:
    """Raise an error if the combination of parameters is invalid"""

    # legacy weights needs classify

    if ctx.obj["legacy_weights"] and not ctx.obj["classify"]:
        logging.error(
            "You have selected --legacy-weights but no classification method. "
            "Please select any --classify method, or remove this parameter"
        )
        raise click.UsageError(
            "You have selected --legacy-weights but no classification method. "
            "Please select any --classify method, or remove this parameter"
        )


def validate_pfam_path(ctx) -> None:
    """Validates whether a BiG-SCAPE db exists when pfam_path is not provided,
    which requires already processed gbk files and hence a DB in output"""

    if ctx.obj["pfam_path"] is None and ctx.obj["db_path"] is None:
        logging.error(
            "Missing option '-p/--pfam-path'."
            "BiG-SCAPE database not provided, a pfam file is "
            "required in order to detect domains."
        )
        raise click.UsageError(
            "Missing option '-p/--pfam-path'."
            "BiG-SCAPE database not provided, a pfam file is "
            "required in order to detect domains."
        )


def validate_domain_include_list(ctx) -> None:
    """Raise an error if both domain includelists are given at the same time"""

    if (
        ctx.obj["domain_includelist_all_path"]
        and ctx.obj["domain_includelist_any_path"]
    ):
        logging.error(
            "You have selected both all and any domain-includelist options. "
            "Please select only one of the two at a time."
        )
        raise click.UsageError(
            "You have selected both all and any domain-includelist options. "
            "Please select only one of the two at a time."
        )


def validate_record_type(ctx, _, record_type) -> Optional[bs_enums.genbank.RECORD_TYPE]:
    """Validates whether a record-type is provided when running classify"""
    valid_types = {mode.value: mode for mode in bs_enums.genbank.RECORD_TYPE}

    if record_type not in valid_types:
        logging.error("Provided --record-type is invalid")
        raise click.UsageError("Provided --record-type in invalid")

    return valid_types[record_type]


def validate_query_record(ctx) -> None:
    """Validates whether a query record number is provided when running query mode
    with a given record type"""

    if (
        ctx.obj["query_record_number"] is None
        and ctx.obj["record_type"] != bs_enums.genbank.RECORD_TYPE.REGION
    ):
        logging.error(
            "Missing option '--query-record-number'."
            "A query record number is required when running query mode with "
            "a record type other than 'region'."
        )
        raise click.UsageError(
            "Missing option '--query-record-number'."
            "A query record number is required when running query mode with "
            "a record type other than 'region'."
        )

    return None
