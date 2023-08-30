"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from pathlib import Path
from typing import List, Optional
from math import floor

# from other modules
from src.genbank.gbk import GBK, SOURCE_TYPE


def load_dataset_folder(
    path: Path,
    source_type: SOURCE_TYPE,
    mode: str = "recursive",
    include_gbk: Optional[List[str]] = None,
    exclude_gbk: Optional[List[str]] = None,
    cds_overlap_cutoff: Optional[float] = None,
    legacy_mode=False,
) -> List[GBK]:
    """Loads all gbk files in a given folder

    Returns empty list if path does not point to a folder or if folder does not contain gbk files
    """

    if not include_gbk:
        include_gbk = ["cluster", "region"]

    if not exclude_gbk:
        exclude_gbk = ["final"]

    if not path.is_dir():
        logging.error("Dataset folder does not point to a directory!")
        raise NotADirectoryError()

    if mode == "recursive":
        files = list(path.glob("**/*.gbk"))

    if mode == "flat":
        files = list(path.glob("*.gbk"))

    # empty folder?
    if len(files) == 0:
        logging.error("Folder does not contain any GBK files!")
        raise FileNotFoundError()

    filtered_files = filter_files(files, include_gbk, exclude_gbk)
    num_files = len(filtered_files)

    logging.info(f"Loading {num_files} GBKs")

    gbk_list = []
    for idx, file in enumerate(filtered_files):
        if num_files > 9 and idx % floor(num_files / 10) == 0:
            logging.info(f"Loaded {idx}/{num_files} GBKs")

        gbk = load_gbk(file, source_type, cds_overlap_cutoff, legacy_mode)

        gbk_list.append(gbk)

    return gbk_list


def filter_files(files: List[Path], include_string: List[str], exclude_str: List[str]):
    """Removes files from input based on include/exclude string conditions

    Args:
        files (List[Path]): list of gbk file paths
        include_string (List[str]): strings that a file must have to be included
        exclude_str (List[str]): strings that a file can't have to be included
    """
    filtered_files = []

    if include_string[0] == "*":
        return files

    for file in files:
        if is_included(file, exclude_str):
            continue
        if is_included(file, include_string):
            filtered_files.append(file)

    return filtered_files


def is_included(path: Path, include_list: List[str]):
    """Returns true if filename includes string from list

    Args:
        path (Path): Path to gbk file
        include_list (List[str]): list of strings to check for presence in filename

    Returns:
        Bool: True if any of strings in list is present in filename
    """
    for string in include_list:
        if string in path.name:
            return True
    return False


def load_gbk(
    path: Path,
    source_type: SOURCE_TYPE,
    cds_overlap_cutoff: Optional[float] = None,
    legacy_mode=False,
) -> GBK:
    """Loads a GBK file. Returns a GBK object

    Args:
        path (Path): path to gbk file
        source_type (SOURCE_TYPE): str, type of gbk file (query, mibig, reference)
        cds_overlap_cutoff (Optional[float]): cds region overlap cutoff to use.
        Defaults to None

    Raises:
        IsADirectoryError: expected file path, got directory instead

    Returns:
        GBK: gbk object
    """

    if not path.is_file():
        logging.error("GBK path does not point to a file!")
        raise IsADirectoryError()

    return GBK.parse(path, source_type, cds_overlap_cutoff, legacy_mode)
