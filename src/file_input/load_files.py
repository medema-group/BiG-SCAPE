"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from pathlib import Path
from typing import List, Optional

# from other modules
from src.genbank.gbk import GBK, CDS, SOURCE_TYPE


def load_dataset_folder(
    path: Path,
    source_type: SOURCE_TYPE,
    mode: str = "recursive",
    include_gbk: Optional[List[str]] = None,
    exclude_gbk: Optional[List[str]] = None,
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

    gbk_list = []
    for file in filtered_files:
        gbk = load_gbk(file, source_type)

        # Filter out CDS for this gbk where CDS coordinates overlap by a certain amount.
        # The longest CDS will be chosen and other CDS will be discarded. This is done
        # in order to remove isoforms of the same gene in organisms that perform
        # alternative splicing
        # TODO: expose percentage to user
        # TODO: do not filter same strand
        gbk_filter_cds_overlap(gbk)

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


def load_gbk(path: Path, source_type: SOURCE_TYPE) -> GBK:
    """Loads a GBK file. Returns a GBK object"""
    # TODO: fix documentation
    if not path.is_file():
        logging.error("GBK path does not point to a file!")
        raise IsADirectoryError()

    return GBK.parse(path, source_type)


def gbk_filter_cds_overlap(gbk: GBK):
    # TODO: document
    # TODO: remove this once the optional problems are gone
    valid_genes = []
    for gene in gbk.genes:
        if gene is None:
            raise ValueError()
        valid_genes.append(gene)

    CDS.filter_overlap(valid_genes)
