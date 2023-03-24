"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from pathlib import Path
from typing import List

# from other modules
from src.genbank.gbk import GBK

# TODO: implement include_gbk/exclude_gbk logics here and account for '*' conflicts


def load_datset_folder(path: Path, mode: str = "recursive") -> List[GBK]:
    """Loads all gbk files in a given folder

    Returns empty list if path does not point to a folder or if folder does not contain gbk files
    """
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

    gbk_list = []
    for file in files:
        gbk = load_gbk(file)

        gbk_list.append(gbk)

    return gbk_list


def load_gbk(path: Path) -> GBK:
    """Loads a GBK file. Returns a GBK object"""
    if not path.is_file():
        logging.error("GBK path does not point to a file!")
        raise IsADirectoryError()

    return GBK.parse(path)
