"""Module to handle loading input GBK data for BiG-SCAPE dereplicate"""

# from python
import logging
from pathlib import Path
import os
import glob
from collections.abc import Iterator, Tuple
import hashlib

# from dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# from other modules
import big_scape.enums as bs_enums
from big_scape.file_input.load_files import filter_files


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
            Path(f) for f in glob.glob(os.path.join(input_dir, "**/*.gbk"), recursive=True)
        ]

    if input_mode == bs_enums.INPUT_MODE.FLAT:
        input_files = list(input_dir.glob("*.gbk"))

    # another redundancy check which should not happen if CLI works correclty
    if len(input_files) == 0:
        logging.error("No GBK files found in the input folder!")
        raise FileNotFoundError()

    filtered_input_files = filter_files(input_files, include_gbk, exclude_gbk)
    num_files = len(filtered_input_files)

    logging.info("Found %s GBK files with valid filenames in the input folder", num_files)

    return filtered_input_files


def parse_gbk_files(input_paths: list[Path]) -> Iterator[Tuple[Path, str, SeqRecord]]:
    """Parse GenBank files and yield their paths, content hashes and SeqIO records

    Args:
        input_paths (list[Path]): input GBK paths

    Yields:
        Iterator[Tuple[Path, str, SeqRecord]]: path, content hash, SeqIO record
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
            logging.warning("%s: GBK contains multiple sequence records, "
                            "only the first one will be used.", path)

        record: SeqRecord = records.pop(0)

        yield (path, hash, record)

