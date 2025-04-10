"""Module to handle loading input GBK data for BiG-SCAPE dereplicate"""

# from python
import logging
from pathlib import Path
import os
import glob
from collections.abc import Iterator
import hashlib
from typing import Optional

# from dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# from other modules
import big_scape.enums as bs_enums
from big_scape.file_input.load_files import filter_files
from big_scape.dereplicating.gbk_component_parsing import get_parser_functions, validate_cds_component
from big_scape.dereplicating.gbk_components.gbk import GBK


def load_input_folder(run: dict) -> list[Path]:
    """Load all valid GBKs files from the input folder and return a list of paths

    Args:
        run (dict): run parameters

    Raises:
        NotADirectoryError: input path is not a directory
        FileNotFoundError: no GBK files found

    Returns:
        list[Path]: valid GBK paths in the input folder
    """

    input_dir = run["input_dir"]
    input_mode = run["input_mode"]
    include_gbk = run["include_gbk"]
    exclude_gbk = run["exclude_gbk"]
    input_files = []

    # redundancy checks, should not happen unless there are issues with the CLI
    if not include_gbk:
        include_gbk = ["cluster", "region"]

    if not exclude_gbk:
        exclude_gbk = ["final"]

    if not input_dir.is_dir():
        logging.error("Dataset folder does not point to a directory!")
        raise NotADirectoryError()

    if input_mode == bs_enums.INPUT_MODE.RECURSIVE:
        input_files = [
            Path(f)
            for f in glob.glob(os.path.join(input_dir, "**/*.gbk"), recursive=True)
        ]

    if input_mode == bs_enums.INPUT_MODE.FLAT:
        input_files = list(input_dir.glob("*.gbk"))

    # another redundancy check which should not happen if CLI works correclty
    if len(input_files) == 0:
        logging.error("No GBK files found in the input folder!")
        raise FileNotFoundError()

    filtered_input_files = filter_files(input_files, include_gbk, exclude_gbk)
    num_files = len(filtered_input_files)

    logging.info(
        "Found %s GBK files with valid filenames in the input folder", num_files
    )

    return filtered_input_files


def parse_gbk_files(
    input_paths: list[Path], source_type: bs_enums.SOURCE_TYPE
) -> Iterator[tuple[Path, str, SeqRecord, bs_enums.SOURCE_TYPE]]:
    """Parse GenBank files and yield their paths, content hashes and SeqIO records

    Args:
        input_paths (list[Path]): input GBK paths
        source_type (bs_enums.SOURCE_TYPE): source type of the GBK

    Yields:
        Iterator[Tuple[Path, str, SeqRecord, bs_enums.SOURCE_TYPE]]: path, content hash, SeqIO record, source type
    """

    for path in input_paths:

        # get unique content hash
        f = open(path, "r")
        content = f.read()
        f.close()
        content = content.encode("utf-8")  # type: ignore
        hash = hashlib.sha256(content).hexdigest()  # type: ignore

        # get SeqIO record
        records = list(SeqIO.parse(path, "genbank"))
        if len(records) > 1:
            logging.warning(
                "%s: GBK contains multiple sequence records, "
                "only the first one will be used.",
                path,
            )

        record: SeqRecord = records.pop(0)

        yield (path, hash, record, source_type)


def gbk_factory(gbk_data: tuple[Path, str, SeqRecord], run: dict) -> Optional[GBK]:
    """Factory function to create a GBK object with all its components

    Args:
        gbk_data (tuple[Path, str, SeqRecord]): GBK file path, hash and SeqIO record
        run (dict): _description_

    Returns:
        GBK: GBK object component
    """

    # get relevant run parameters
    run_mode = run["mode"]

    path, hash, seqIO_record, source_type = gbk_data

    # we dont need to store the entire nt_sequence, just the length
    # so we can filter out the GBKs that are too long/short later on
    nt_length = len(seqIO_record.seq)

    # get antiSMASH version
    as_version = get_as_version(seqIO_record)

    # create GBK object
    gbk: GBK = GBK.create(path, hash, nt_length, as_version, source_type)

    gbk = parse_seqIO(gbk, seqIO_record, run_mode)

    if not validate_cds_component(gbk):
        logging.error("This GBK file does not contain any CDS features, %s", path)
        return None

    # TODO: hybrid collapsing for protoclusters

    # TODO: if run mode is not derep and force-gbk = True, create fake region if region not there 
    #       if region not there and force-gbk = False, raise error

    # TODO validate that all the parenting is there properly (follow top down hierarchy)

    return gbk


def parse_seqIO(gbk: GBK, seqIO_record: SeqRecord, run_mode) -> GBK:
    """Parse SeqIO record and add components to GBK object

    Args:
        gbk (GBK): gbk object
        seqIO_record (SeqRecord): SeqIO record

    Returns:
        GBK: gbk object with components added
    """

    # get relevant function dict for feature parsing

    parser_functions = get_parser_functions(run_mode)

    # parse SeqIO record and add components to GBK object

    for feature in seqIO_record.features:

        feature_type = feature.type

        # skip features we dont ever want to parse
        if feature_type not in set(feature.value for feature in bs_enums.FEATURE_TYPE):
            continue

        feature_name = bs_enums.FEATURE_TYPE(feature_type)

        # skip features we dont care about
        if feature_name not in parser_functions:
            continue

        # get parser function for feature
        parser_function = parser_functions[feature_name]

        component = parser_function(feature, seqIO_record, gbk)

        if not component:
            continue

        # remember that cluster type features are returned as region type components

        gbk.components.setdefault(component.name, []).append(component)

    return gbk


def get_as_version(gbk_seq_record: SeqRecord) -> str:
    """Get AS version from GBK record

    Args:
        gbk_seq_record (SeqRecord): gbk seqrecord

    Returns:
        str: antismash version
    """

    try:
        as_version = gbk_seq_record.annotations["structured_comment"][
            "antiSMASH-Data"
        ]["Version"]
    except KeyError:
        # assume AS version 4
        as_version = "4"

    return as_version