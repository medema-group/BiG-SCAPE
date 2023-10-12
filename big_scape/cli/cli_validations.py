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

    param_dict["label"] = f"{timestamp}"

    param_dict["start_time"] = start_time

    return None


# meta parameter validations


def validate_profiling(ctx, param, profiling) -> bool:
    """Checks whether multithreading is possible, and therefore whether profiling can happen"""

    if profiling and platform.system() == "Darwin":
        logging.warning("Profiling is not supported on MacOS, please use Linux")

    return False


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


# output parameter validations


def validate_output_paths(ctx) -> None:
    """Sets the output paths to default output_dir if not provided"""

    timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())

    if ctx.obj["db_path"] is None:
        db_path = ctx.obj["output_dir"] / Path("data_sqlite.db")
        ctx.obj["db_path"] = db_path

    if ctx.obj["log_path"] is None:
        log_filename = timestamp + ".log"
        log_path = ctx.obj["output_dir"] / Path(log_filename)
        ctx.obj["log_path"] = log_path

    if ctx.obj["profile_path"] is None:
        profile_filename = timestamp + ".profile"
        profile_path = ctx.obj["output_dir"] / Path(profile_filename)
        ctx.obj["profile_path"] = profile_path

    return None


# comparison validations


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

    return gcf_cutoffs


def validate_filter_gbk(ctx, param, filter_str) -> list[str]:
    """Validates and formats the filter string and returns a list of strings"""

    return filter_str.split(",")


# hmmer parameters


def validate_includelist(ctx, param, domain_includelist_path) -> None:
    """Validate the path to the domain include list and return a list of domain
    accession strings contained within this file

    Returns:
        list[str]: A list of domain accessions to include
    """

    # only validate if set
    if domain_includelist_path is None:
        return None

    if not domain_includelist_path.exists():
        logging.error("domain_includelist file does not exist!")
        raise InvalidArgumentError("--domain_includelist", domain_includelist_path)

    with domain_includelist_path.open(encoding="utf-8") as domain_includelist_file:
        lines = domain_includelist_file.readlines()

        lines = [line.strip() for line in lines]

        # expect Pfam accessions, i.e. PF00001 or PF00001.10
        lines_valid = map(
            lambda string: string.startswith("PF") and len(string) in range(7, 11),
            lines,
        )

        if not all(lines_valid):
            logging.error(
                "Invalid Pfam accession(s) found in file %s", domain_includelist_path
            )
            raise click.BadParameter(
                "Invalid Pfam accession(s) found in file %s", domain_includelist_path
            )

        ctx.params["domain_includelist"] = lines


# workflow validations


def validate_binning_workflow(ctx) -> None:
    """Raise an error if the combination of parameters in this object means no work
    will be done
    """

    if (
        ctx.obj["no_mix"] is True
        and ctx.obj["legacy_classify"] is False
        and ctx.obj["classify"] is False
    ):
        logging.error(
            "The combination of arguments you have selected for binning means no work will "
            "be done. Please remove either --no_mix, or add --legacy_classify/--classify"
            " in order to enable comparisons"
        )
        raise click.UsageError(
            "The combination of arguments you have selected for binning means no work will "
            "be done. Please remove either --no_mix, or add --legacy_classify/--classify"
            " in order to enable comparisons"
        )


def validate_skip_hmmscan(ctx) -> None:
    """Validates whether a BiG-SCAPE db exists when running skip_hmm, which
    requires already processed gbk files and hence a DB in output"""

    if ctx.obj["skip_hmmscan"] and ctx.obj["db_path"] is None:
        logging.error(
            "Missing option '--db_path'."
            "BiG-SCAPE database has not been given, skip_hmmscan requires "
            "a DB of already processed gbk files."
        )
        raise click.UsageError(
            "Missing option '--db_path'."
            "BiG-SCAPE database has not been given, skip_hmmscan requires "
            "a DB of already processed gbk files."
        )


def validate_pfam_path(ctx) -> None:
    """Validates whether a BiG-SCAPE db exists when pfam_path is not provided,
    which requires already processed gbk files and hence a DB in output"""

    if ctx.obj["pfam_path"] is None and ctx.obj["db_path"] is None:
        logging.error(
            "Missing option '-p/--pfam_path'."
            "BiG-SCAPE database not provided, a pfam file is "
            "required in order to detect domains."
        )
        raise click.UsageError(
            "Missing option '-p/--pfam_path'."
            "BiG-SCAPE database not provided, a pfam file is "
            "required in order to detect domains."
        )