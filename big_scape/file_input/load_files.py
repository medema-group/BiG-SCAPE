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
from big_scape.genbank.gbk import GBK

import big_scape.parameters as bs_param
import big_scape.enums as bs_enums
import big_scape.genbank as bs_gbk
import big_scape.data as bs_data
import big_scape.hmm as bs_hmm


def get_mibig(mibig_version: str, bigscape_dir: Path):
    # pragma: no cover
    """A function to download a given version of MIBiG (if need be), and return its location
        pragma: no cover -> we dont wish to test downloading and decompressing of files, due to
        time and resource constraints

    Args:
        mibig_version (str): version
        bigscape_dir (Path): main BiG-SCAPE directory

    Returns:
        Path: path to MIBiG database (antismash processed gbks)
    """

    mibig_url_base = "https://dl.secondarymetabolites.org/mibig/"
    mibig_url = mibig_url_base + f"mibig_antismash_{mibig_version}_gbk.tar.bz2"
    # TODO: this only works for 3.1, update to proper link once Kai makes it available
    # https://dl.secondarymetabolites.org/mibig/mibig_antismash_3.1_gbk.tar.bz2
    # https://dl.secondarymetabolites.org/mibig/mibig_antismash_3.1_json.tar.bz2

    mibig_dir = Path(os.path.join(bigscape_dir, "MIBiG"))
    mibig_version_dir = Path(
        os.path.join(mibig_dir, f"mibig_antismash_{mibig_version}_gbk")
    )

    if not os.path.exists(mibig_dir):
        os.makedirs(mibig_dir)

    contents_dir = os.listdir(mibig_dir)
    if len(contents_dir) > 0 and any(mibig_version in dir for dir in contents_dir):
        # we assume that if a folder is here, that it is uncompressed and ready to use
        return mibig_version_dir

    mibig_dir_compressed = Path(f"{mibig_version_dir}.tar.bz2")
    download_dataset(mibig_url, mibig_dir, mibig_dir_compressed)
    return mibig_version_dir


def download_dataset(url: str, path: Path, path_compressed: Path) -> None:
    # pragma: no cover
    """A function to download and decompress a dataset from an online repository
        pragma: no cover -> we dont wish to test downloading and decompressing of files, due to
        time and resource constraints

    Args:
        url (str): url to download
        path (Path): path to decompressed file/folder
        path_compressed (Path): path to compressed file/folder

    Raises:
        ValueError: could not download/request file

    Returns:
        None
    """

    # download
    response = requests.get(url, stream=True, timeout=10)

    if response.status_code != 200:
        raise ValueError("Could not download MIBiG file")

    with open(path_compressed, "wb") as file:
        file.write(response.raw.read())

    # extract contents
    with tarfile.open(path_compressed) as file:
        file.extractall(path)

    os.remove(path_compressed)


def load_dataset_folder(
    path: Path,
    source_type: bs_enums.SOURCE_TYPE,
    mode: str = "recursive",
    include_gbk: Optional[List[str]] = None,
    exclude_gbk: Optional[List[str]] = None,
    cds_overlap_cutoff: Optional[float] = None,
    legacy_mode=False,
) -> List[GBK]:
    """Loads all gbk files in a given folder

    Returns empty list if path does not point to a folder or if folder does not contain
    gbk files
    """

    if source_type in (bs_enums.SOURCE_TYPE.QUERY, bs_enums.SOURCE_TYPE.REFERENCE):
        if not include_gbk:
            include_gbk = ["cluster", "region"]

        if not exclude_gbk:
            exclude_gbk = ["final"]

    if source_type == bs_enums.SOURCE_TYPE.MIBIG:
        include_gbk = ["BGC"]

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

    logging.info("Loading %d GBKs", num_files)

    gbk_list = []
    for idx, file in enumerate(filtered_files):
        if num_files > 9 and idx % floor(num_files / 10) == 0:
            logging.info("Loaded %d/%d GBKs", idx, num_files)

        gbk = load_gbk(file, source_type, cds_overlap_cutoff, legacy_mode)

        gbk_list.append(gbk)

    return gbk_list


def filter_files(
    files: List[Path],
    include_string: Optional[List[str]] = None,
    exclude_str: Optional[List[str]] = None,
):
    """Removes files from input based on include/exclude string conditions

    Args:
        files (List[Path]): list of gbk file paths
        include_string (List[str]): strings that a file must have to be included
        exclude_str (List[str]): strings that a file can't have to be included
    """
    filtered_files = []

    if include_string and include_string[0] == "*":
        return files

    for file in files:
        if exclude_str and is_included(file, exclude_str):
            continue
        if include_string and is_included(file, include_string):
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
    source_type: bs_enums.SOURCE_TYPE,
    cds_overlap_cutoff: Optional[float] = None,
    legacy_mode=False,
) -> GBK:
    """Loads a GBK file. Returns a GBK object

    Args:
        path (Path): path to gbk file
        source_type (bs_enums.SOURCE_TYPE): str, type of gbk file (query, mibig,
            reference)
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


def load_gbks(run: bs_param.RunParameters, bigscape_dir: Path) -> list[GBK]:
    """Load all gbks from the input folder and return a list of GBK objects

    Args:
        run (bs_param.RunParameters): run parameters
        bigscape_dir (Path): path to BiG-SCAPE directory

    Returns:
        list[GBK]: list of GBK objects
    """

    # counterintuitive, but we need to load the gbks first to see if there are any
    # differences with the data in the database

    input_gbks = []

    if run.binning.query_bgc_path:
        query_bgc_gbk = GBK(run.binning.query_bgc_path, bs_enums.SOURCE_TYPE.QUERY)
        input_gbks.append(query_bgc_gbk)

        query_bgc_stem = run.binning.query_bgc_path.stem

        # add the query bgc to the exclude list
        exclude_gbk = run.input.exclude_gbk + [query_bgc_stem]

        gbks = load_dataset_folder(
            run.input.input_dir,
            bs_enums.SOURCE_TYPE.REFERENCE,
            run.input.input_mode,
            run.input.include_gbk,
            exclude_gbk,
            run.input.cds_overlap_cutoff,
        )
        input_gbks.extend(gbks)

    else:
        gbks = load_dataset_folder(
            run.input.input_dir,
            bs_enums.SOURCE_TYPE.QUERY,
            run.input.input_mode,
            run.input.include_gbk,
            run.input.exclude_gbk,
            run.input.cds_overlap_cutoff,
        )
        input_gbks.extend(gbks)

    # get reference if either MIBiG version or user-made reference dir passed
    if run.input.mibig_version:
        mibig_version_dir = get_mibig(run.input.mibig_version, bigscape_dir)
        mibig_gbks = load_dataset_folder(mibig_version_dir, bs_enums.SOURCE_TYPE.MIBIG)
        input_gbks.extend(mibig_gbks)

    if run.input.reference_dir:
        reference_gbks = load_dataset_folder(
            run.input.reference_dir, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(reference_gbks)

    # find the minimum task set for these gbks
    # if there is no database, create a new one and load in all the input stuff
    if not run.output.db_path.exists():
        bs_data.DB.create_in_mem()

        all_cds: list[bs_gbk.CDS] = []
        for gbk in input_gbks:
            gbk.save_all()
            all_cds.extend(gbk.genes)

        return input_gbks

    bs_data.DB.load_from_disk(run.output.db_path)
    task_state = bs_data.find_minimum_task(input_gbks)

    # if we are are not on the load_gbks task, we have all the data we need
    if task_state != bs_enums.TASK.LOAD_GBKS:
        logging.info("Loading existing run from disk...")

        gbks = GBK.load_all()

        for gbk in gbks:
            bs_hmm.HSP.load_all(gbk.genes)

        return gbks

    # if we end up here, we are in some halfway state and need to load in the new data
    logging.info("Loading existing run from disk and adding new data...")
    missing_gbks = bs_data.get_missing_gbks(input_gbks)
    logging.info("Found %d missing gbks", len(missing_gbks))

    for gbk in missing_gbks:
        gbk.save_all()

    # still return the full set
    return input_gbks
