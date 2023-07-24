"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from pathlib import Path
from typing import List, Optional
from math import floor
import os
import tarfile

# from dependencies
import requests  # type: ignore

# from other modules
from src.genbank.gbk import GBK, SOURCE_TYPE


def get_mibig(mibig_version: str, reference_dir: Path):
    """download mibig"""

    # if there is a given valid path, use (no version, given path)
    if reference_dir and os.path.exists(reference_dir):
        return reference_dir

    # TODO: update to proper link once Kai makes it available
    mibig_url = (
        f"https://dl.secondarymetabolites.org/mibig/mibig_json_{mibig_version}.tar.gz"
    )

    # download to given path (given version, given path but empty)
    if mibig_version and reference_dir and not os.path.exists(reference_dir):
        reference_dir_compressed = Path(f"{reference_dir}.tar.gz")
        download_dataset(mibig_url, reference_dir, reference_dir_compressed)
        return reference_dir

    default_parent_dir = Path(os.path.dirname(os.path.abspath(__file__)), "reference")
    default_reference_dir = Path(default_parent_dir, f"mibig_antismash_{mibig_version}")
    default_reference_dir_compressed = Path(str(reference_dir) + ".tar.gz")

    # download to default path (yes version, no path - default full)
    if mibig_version and os.path.exists(default_reference_dir):
        return default_reference_dir

    # download to default path (yes version, no path - default empty)
    if mibig_version and not os.path.exists(default_reference_dir):
        download_dataset(mibig_url, reference_dir, default_reference_dir_compressed)
        return default_reference_dir

    return reference_dir


def get_pfam(pfam_version: Optional[str], pfam_path: Optional[Path]):
    """Download pfam file"""

    # defaults
    def_pfam_dir = Path(os.path.dirname(os.path.abspath(__file__)), "pfam")
    def_pfam_path = Path(def_pfam_dir, "Pfam-A.hmm")
    def_pfam_path_compressed = Path(def_pfam_dir, "Pfam-A.hmm.gz")
    current_version_url = (
        "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    )

    # values to set
    if pfam_version:
        version_url = f"https://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam{pfam_version}/Pfam-A.hmm.gz"

    if pfam_path:
        pfam_path_compressed = Path(f"{pfam_path}.gz")

    # if there is a path given

    if pfam_path:
        # if path already filled, use
        if pfam_path.exists():
            return pfam_path

        if not pfam_path.exists():
            if pfam_version:
                download_dataset(version_url, pfam_path, pfam_path_compressed)
                return pfam_path

            else:
                download_dataset(current_version_url, pfam_path, pfam_path_compressed)
                return pfam_path

    # if there is no path given

    if not pfam_path:
        # if default path already filled, use it
        if def_pfam_path.exists():
            return def_pfam_path

        # if its the first time and have to make the parent dir
        if not def_pfam_dir.exists():
            os.makedirs(def_pfam_dir)

        # download given or current version
        if pfam_version:
            download_dataset(version_url, def_pfam_path, def_pfam_path_compressed)
        else:
            download_dataset(
                current_version_url, def_pfam_path, def_pfam_path_compressed
            )
        return def_pfam_path


def download_dataset(url: str, path: Path, path_compressed: Path) -> None:
    """Download a dataset"""

    response = requests.get(url, stream=True)
    if response.status_code != 200:
        raise ValueError("Could not download MIBiG file")
    with open(path_compressed, "wb") as f:
        f.write(response.raw.read())

    # extract contents
    file = tarfile.open(path_compressed)
    file.extractall(path)
    file.close()
    os.remove(path_compressed)

    return None


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
