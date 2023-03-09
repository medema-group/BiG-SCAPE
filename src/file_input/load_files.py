"""Module to handle automatically loading files (iterating directories, etc)"""

import logging

from pathlib import Path

from src.genbank.gbk import GBK


def load_datset_folder(path: Path) -> list[GBK]:
    """Loads all gbk files in a given folder

    Returns empty list if path does not point to a folder or if folder does not contain gbk files
    """
    if not path.is_dir():
        logging.error("Dataset folder does not point to a directory!")
        raise NotADirectoryError()

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

    return GBK.parse_gbk(path)
