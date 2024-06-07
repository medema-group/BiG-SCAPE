"""Module to handle automatically loading files (iterating directories, etc)"""

# from python
import logging
from pathlib import Path
from typing import List, Optional
import os
import tarfile
import multiprocessing

# from dependencies
import requests  # type: ignore

# from other modules
from big_scape.genbank.gbk import GBK
from big_scape.cli.config import BigscapeConfig
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
        logging.info("Creating MIBiG %s directory", mibig_version)
        os.makedirs(mibig_dir)

    contents_dir = os.listdir(mibig_dir)
    if len(contents_dir) > 0 and any(mibig_version in dir for dir in contents_dir):
        logging.info("MIBiG version %s already downloaded", mibig_version)
        # we assume that if a folder is here, that it is uncompressed and ready to use
        return mibig_version_dir

    logging.info("Downloading MIBiG version %s", mibig_version)
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
        # TODO: deal with deprecation
        file.extractall(path)

    os.remove(path_compressed)


def load_dataset_folder(
    input_dir: Path, run: dict, source_type: bs_enums.SOURCE_TYPE
) -> List[GBK]:
    """Loads all gbk files in a given folder

    Returns empty list if path does not point to a folder or if folder does not contain
    gbk files
    """

    path = input_dir
    source_type = source_type
    mode = run["input_mode"]
    include_gbk = run["include_gbk"]
    exclude_gbk = run["exclude_gbk"]
    cds_overlap_cutoff = BigscapeConfig.CDS_OVERLAP_CUTOFF
    cores = run["cores"]

    if source_type in (bs_enums.SOURCE_TYPE.QUERY, bs_enums.SOURCE_TYPE.REFERENCE):
        # redundant check, this should never happen if the cli is working properly
        if not include_gbk:
            include_gbk = ["cluster", "region"]

        if not exclude_gbk:
            exclude_gbk = ["final"]

    if source_type == bs_enums.SOURCE_TYPE.MIBIG:
        include_gbk = ["BGC"]

    if not path.is_dir():
        logging.error("Dataset folder does not point to a directory!")
        raise NotADirectoryError()

    if mode == bs_enums.INPUT_MODE.RECURSIVE:
        files = list(path.glob("**/*.gbk"))

    if mode == bs_enums.INPUT_MODE.FLAT:
        files = list(path.glob("*.gbk"))

    # redundant check, this should never happen if the cli is working properly
    if len(files) == 0:
        logging.error("Folder does not contain any GBK files!")
        raise FileNotFoundError()

    filtered_files = filter_files(files, include_gbk, exclude_gbk)
    num_files = len(filtered_files)

    logging.info("Loading %d GBKs", num_files)

    if cores is None:
        cores = multiprocessing.cpu_count()

    gbk_list = []
    pool = multiprocessing.Pool(cores)

    # TODO: would really like to see a loading bar here
    gbk_list = pool.starmap(
        load_gbk,
        [(file, source_type, run, cds_overlap_cutoff) for file in filtered_files],
    )

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


def remove_duplicate_gbk(gbks: list[GBK]) -> list[GBK]:
    """Remove any duplicate GBKs from given input based on content hash

    Args:
        gbks (list[GBK]): list of input GBKs

    Returns:
        list[GBK]: list of unique input GBKs
    """
    unique_gbks = set()
    for gbk in gbks:
        if gbk in unique_gbks:
            logging.info("Skipping duplicate %s", gbk)
            continue
        unique_gbks.add(gbk)

    return list(unique_gbks)


def bgc_length_contraint(gbks: list[GBK]) -> list[GBK]:
    """Remove any GBKs that are not between the minimum ans maximum allowed length

    Args:
        gbks (list[GBK]): list of input GBKs

    Returns:
        list[GBK]: list of filtered input GBKs
    """
    orig_size = len(gbks)
    filtered_gbks = [
        gbk
        for gbk in gbks
        if BigscapeConfig.MIN_BGC_LENGTH
        < len(gbk.nt_seq)
        < BigscapeConfig.MAX_BGC_LENGTH
    ]
    filtered_size = len(filtered_gbks)

    if filtered_size == 0:
        logging.error(
            "No GBKs remain after length constraint of min %s to max %s bp!",
            BigscapeConfig.MIN_BGC_LENGTH,
            BigscapeConfig.MAX_BGC_LENGTH,
        )
        raise RuntimeError()

    if orig_size != filtered_size:
        logging.info(
            "%s GBKs remain after length constraint of min %s to max %s bp",
            filtered_size,
            BigscapeConfig.MIN_BGC_LENGTH,
            BigscapeConfig.MAX_BGC_LENGTH,
        )
    return filtered_gbks


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
    run: dict,
    cds_overlap_cutoff: Optional[float] = None,
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

    return GBK.parse(path, source_type, run, cds_overlap_cutoff)


def load_gbks(run: dict, bigscape_dir: Path) -> list[GBK]:
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

    if run["query_bgc_path"]:
        query_bgc_gbk = GBK.parse(
            run["query_bgc_path"],
            bs_enums.SOURCE_TYPE.QUERY,
            run,
            run["cds_overlap_cutoff"],
        )
        logging.info("Loaded query GBK %s", query_bgc_gbk)

        input_gbks.append(query_bgc_gbk)

        query_bgc_stem = run["query_bgc_path"].stem

        # add the query bgc to the exclude list
        run["exclude_gbk"].append(query_bgc_stem)

        gbks = load_dataset_folder(
            run["input_dir"], run, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(gbks)

    else:
        gbks = load_dataset_folder(run["input_dir"], run, bs_enums.SOURCE_TYPE.QUERY)
        input_gbks.extend(gbks)

    # get reference if either MIBiG version or user-made reference dir passed
    if run["mibig_version"]:
        mibig_gbks = load_dataset_folder(
            run["mibig_dir"], run, bs_enums.SOURCE_TYPE.MIBIG
        )
        input_gbks.extend(mibig_gbks)

    if run["reference_dir"]:
        reference_gbks = load_dataset_folder(
            run["reference_dir"], run, bs_enums.SOURCE_TYPE.REFERENCE
        )
        input_gbks.extend(reference_gbks)

    # remove eventual duplicates from input gbks
    input_gbks = remove_duplicate_gbk(input_gbks)

    # apply minimum and maximum bgc length constraint
    input_gbks = bgc_length_contraint(input_gbks)

    # find the minimum task set for these gbks
    task_state = bs_data.find_minimum_task(input_gbks)

    source_dict = {gbk.hash: gbk.source_type for gbk in input_gbks}

    # if we are are not on the save_gbks task, we have all the data we need in the database
    # and can just load it all into the correct python objects
    if task_state != bs_enums.TASK.SAVE_GBKS:
        # here we dont save anything to DB, data goes DB -> python objects
        logging.info("Loading existing run from disk...")

        input_gbks_from_db = GBK.load_many(input_gbks)
        for gbk in input_gbks_from_db:
            gbk.source_type = source_dict[gbk.hash]
            bs_hmm.HSP.load_all(gbk.genes)

        return input_gbks_from_db

    # if we end up here, we are in some halfway state and need to load in the new data
    logging.info("Loading existing run from disk and adding new data...")
    missing_gbks = bs_data.get_missing_gbks(input_gbks)
    logging.info("Found %d missing gbks", len(missing_gbks))

    for gbk in missing_gbks:
        gbk.save_all()

    # now we have all new data in the database, we can load it all in to the correct
    # python objects
    input_gbks_from_db = GBK.load_many(input_gbks)
    for gbk in input_gbks_from_db:
        gbk.source_type = source_dict[gbk.hash]
        bs_hmm.HSP.load_all(gbk.genes)

    return input_gbks_from_db


def get_all_bgc_records(run: dict, gbks: List[GBK]) -> List[bs_gbk.BGCRecord]:
    """Get all BGC records from the working list of GBKs

    Args:
        gbks (list[GBK]): list of GBK objects
        run (dict): run parameters

    Returns:
        list[bs_gbk.BGCRecord]: list of BGC records
    """
    all_bgc_records: list[bs_gbk.BGCRecord] = []
    for gbk in gbks:
        if gbk.region is not None:
            gbk_records = bs_gbk.bgc_record.get_sub_records(
                gbk.region, run["record_type"]
            )
            all_bgc_records.extend(gbk_records)

    return all_bgc_records


def get_all_bgc_records_query(
    run: dict, gbks: List[GBK]
) -> tuple[List[bs_gbk.BGCRecord], bs_gbk.BGCRecord]:
    """Get all BGC records from the working list of GBKs

    Args:
        gbks (list[GBK]): list of GBK objects
        run (dict): run parameters

    Returns:
        list[bs_gbk.BGCRecord]: list of BGC records
    """
    all_bgc_records: list[bs_gbk.BGCRecord] = []
    for gbk in gbks:
        if gbk.region is not None:
            gbk_records = bs_gbk.bgc_record.get_sub_records(
                gbk.region, run["record_type"]
            )
            if gbk.source_type == bs_enums.SOURCE_TYPE.QUERY:
                query_record_type = run["record_type"]

                query_record_type = run["record_type"]
                query_record_number = run["query_record_number"]

                query_sub_records = bs_gbk.bgc_record.get_sub_records(
                    gbk.region, query_record_type
                )

                if query_record_type == bs_enums.RECORD_TYPE.REGION:
                    query_record = query_sub_records[0]

                else:
                    query_record = [
                        record
                        for record in query_sub_records
                        if record.number == query_record_number
                    ][0]

                all_bgc_records.append(query_record)

            else:
                all_bgc_records.extend(gbk_records)

    return all_bgc_records, query_record
