"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from itertools import combinations
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

        filter_cds_overlap(gbk)

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


def filter_cds_overlap(gbk: GBK):
    # TODO: document
    # TODO: remove this once the optional problems are gone
    valid_genes = []
    for gene in gbk.genes:
        if gene is None:
            raise ValueError()
        valid_genes.append(gene)

    del_list = set()
    # find all combinations of cds to check for overlap
    cds_a: CDS
    cds_b: CDS
    for cds_a, cds_b in combinations(valid_genes, 2):
        a_start = int(cds_a.nt_start)
        b_start = int(cds_b.nt_start)

        a_stop = int(cds_a.nt_stop)
        b_stop = int(cds_b.nt_stop)

        # these are different from calculating end - start / 3
        # or something similar.
        # the number that results from this is what is used in the
        # original implementation of this functionality, so the same is
        # done here
        a_len = len(cds_a.aa_seq)
        b_len = len(cds_b.aa_seq)

        # copy from original
        if b_stop <= a_start or b_start >= a_stop:
            pass
        else:
            # calculate overlap
            if a_start > b_start:
                overlap_start = a_start
            else:
                overlap_start = b_start

            if a_stop < b_stop:
                overlap_end = a_stop
            else:
                overlap_end = b_stop

            overlap_length = overlap_end - overlap_start

            # allow the overlap to be as large as 10% of the
            # shortest CDS. Overlap length is in nucleotides
            # here, whereas a_len, b_len are protein
            # sequence lengths
            if overlap_length / 3 > 0.1 * min(a_len, b_len):
                if a_len > b_len:
                    del_list.add(cds_b)
                else:
                    del_list.add(cds_a)
    # remove any entries that need to be removed
    for cds in del_list:
        gbk.genes.remove(cds)
