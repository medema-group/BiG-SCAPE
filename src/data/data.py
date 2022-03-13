"""Module containing general helper functions to load input data into the database
"""

import logging
import os
import re
import glob

from os import path
from multiprocessing import Pool

from src.big_scape import Run
from src.data.hmm import load_hmms
from .database import Database
from .bgc import BGC

def parse_input_gbk(arguments: tuple):
    """Unwraps argument tuple and parse GBK"""
    orig_gbk_path, file_path = arguments
    return (file_path, BGC.parse_gbk(file_path, orig_gbk_path=orig_gbk_path))

def insert_dataset(database, dataset_name, dataset_meta):
    """This function inserts dataset source information into the database
    If the datasets already exist, this will return the stored bgc ids per dataset
    in a set. Otherwise, it wil return an empty set"""
    # This only store bgcs from gbk files
    # exist within the folder
    # e.g. if a dataset has been parsed before,
    # then a folder is deleted, the previous BGCs
    # won't be included here
    bgc_ids = set()

    # check if dataset exists in database
    queried_dataset = database.select(
        "dataset",
        "WHERE name=?",
        parameters=(dataset_name,),
        props=["id"])
    assert len(queried_dataset) <= 1

    if len(queried_dataset) > 0:
        # dataset exists, use entries from db
        # note:
        # in the future, we would like to match
        # between what's in the database
        # and what's in the folder, and adjust
        # the included BGCs accordingly
        # i.e. if there's a new file to parse
        dataset_id = queried_dataset[0]["id"]
        bgc_ids = [
            row["id"] for row in database.select(
                "bgc",
                "WHERE dataset_id=?",
                parameters=(dataset_id,),
                props=["id"])]
        logging.info("Found %d BGCs in the database.", len(bgc_ids))
    else:  # create a new dataset entry
        dataset_id = database.insert(
            "dataset",
            {
                "name": dataset_name,
                "orig_folder": dataset_meta["path"],
                "description": dataset_meta["desc"]
            }
        )

    return dataset_id, bgc_ids

def list_gbk_files(data_path):
    """Returns a list of full paths to gbk files under a certain directory
    This lists files recursively (also in directories of directories)
    """
    files = []
    # define eligible regexes for clustergbks
    eligible_regexes = [re.compile(rgx) for rgx in [
        "^BGC[0-9]{7}\.?[0-9]?$",  # MIBiG
        "^.+\\.cluster[0-9]+$",  # antiSMASH4 clustergbks
        "^.+\\.region[0-9]+$",  # antiSMASH5 clustergbks
    ]]

    # fetch gbk files
    for gbk_full_path in glob.iglob(path.join(
            data_path,
            "**/*.gbk"),
            recursive=True):

        # first check: only take antiSMASH/MIBiG-derived GBKs
        gbk_name = path.splitext(path.basename(gbk_full_path))[0]
        eligible_file = False
        for eligible_pattern in eligible_regexes:
            if eligible_pattern.match(gbk_name):
                eligible_file = True
                break

        # second check: see if already exists in db
        if eligible_file:
            files.append(gbk_full_path)

    return files

def check_gbk_exists(run, database: Database, dataset_id: int, gbk_path: str):
    """Returns the presence of a given gbk file in the database"""
    gbk_name = path.splitext(gbk_path)[0]
    names = database.select("bgc",
                            f"where dataset_id = {dataset_id}\
                              and name = \"{gbk_name}\"")
    return len(names) > 0

def insert_dataset_gbks(run, database: Database, dataset_id, dataset_name, dataset_meta, bgc_ids):
    """Performs the insertion of GBK information into the database"""
    new_bgcs_count = 0
    files_to_process = []
    count_gbk_exists = 0


    for gbk_full_path in list_gbk_files(dataset_meta["path"]):
        gbk_path = path.join(path.dirname(gbk_full_path).split("/")[-1], path.basename(gbk_full_path))
        # note: this check is turned off for now,
        # see above note

        bgc_ids = set()
        if check_gbk_exists(run, database, dataset_id, gbk_path):
            count_gbk_exists += 1
            bgc_ids.update(bgc_ids)
        else:
            files_to_process.append((gbk_path, gbk_full_path))

    if len(files_to_process) == 0:
        logging.info("Found no new GBK files.")
        return

    if count_gbk_exists > 0:
        logging.info("Found %d existing GBKs...", count_gbk_exists)

    # parse and insert new GBKs #
    logging.info("Parsing and inserting %d new GBKs...", len(files_to_process))
    mp_pool = Pool(run.options.cores)
    pool_results = mp_pool.map(parse_input_gbk, files_to_process)
    for file_path, bgcs in pool_results:
        for bgc in bgcs:
            bgc.save(dataset_id, database)
            new_bgcs_count += 1
            bgc_ids.add(bgc.id)
            database.insert(
                "bgc_status",
                {
                    "bgc_id": bgc.id,
                    "status": 1
                }, True
            )
    database.commit_inserts()
    logging.info("Inserted %d new BGCs.", new_bgcs_count)
    return bgc_ids

def create_bgc_status(db: Database, bgc_ids):
    """Create rows for the bgc_status table for each id passed in bgc_ids
    """
    # dictionary of {bgc_id: 1}
    bgc_ids_dict = {bgc_id: 1 for bgc_id in bgc_ids}
    db.insert("bgc_status", bgc_ids_dict)

def initialize_db(run: Run, database: Database):
    """Fills the database with input data"""
    logging.debug("Initializing database")

    datasets = dict()
    if run.mibig.use_mibig:
        datasets["mibig"] = {
            "path": run.mibig.gbk_path,
            "desc": "Mibig dataset"
        }

    datasets["input"] = {
        "path": run.directories.input,
        "desc": "Input files"
    }

    for dataset_name, dataset_meta in datasets.items():
        dataset_id, bgc_ids = insert_dataset(database, dataset_name, dataset_meta)
        bgc_ids = insert_dataset_gbks(run, database, dataset_id, dataset_name, dataset_meta, bgc_ids)

    load_hmms(run, database)




def get_cluster_id_list(database):
    """Returns a list of all bgc ids in the database"""
    bgc_ids = [
        row["id"] for row in database.select(
            "bgc",
            "",
            props=["id"])]
    return bgc_ids


def get_cluster_name_list(database):
    """Returns a list of all bgc ids in the database"""
    bgc_ids = [
        row["name"] for row in database.select(
            "bgc",
            "",
            props=["name"])]
    return bgc_ids



def get_cluster_gbk_dict(run: Run, database: Database):
    """Gets a list of source gbk file paths for each bgc in the database"""
    # mibig paths
    gbk_dict = dict()
    if run.mibig.use_mibig:
        rows = database.select(
            "bgc",
            "where bgc.dataset_id = (select id from dataset where name = \"mibig\")",
            props=["name", "orig_filename"])
        for row in rows:
            gbk_dict[row["name"]] = os.path.join(run.directories.mibig, row["orig_filename"])

    # input paths
    rows = database.select(
        "bgc",
        "where bgc.dataset_id = (select id from dataset where name = \"input\")",
        props=["name", "orig_folder", "orig_filename"])
    for row in rows:
        gbk_dict[row["name"]] = os.path.join(run.directories.input, row["orig_folder"], row["orig_filename"])

    return gbk_dict

def gen_bgc_info_for_svg(database: Database):
    """Generates a dictionary of BGC info objects that are needed for the SVG image generation"""
    bgc_info = dict()
    rows = database.select(
        "cds",
        "join bgc on cds.bgc_id = bgc.id group by bgc_id",
        props=["bgc_id", "bgc.name as bgc_name", "max(nt_end - nt_start) as max_width", "count(cds.id) as records"])
    for row in rows:
        bgc_info[row["bgc_name"]] = {
            "max_width": row["max_width"],
            "records": row["records"]
        }
    return bgc_info

def gen_bgc_info_for_fetch_genome(database: Database):
    """Generates a dictionary of BGC info objects that are needed for the SVG image generation"""
    bgc_info = dict()
    rows = database.select(
        "bgc",
        props=["bgc.name as bgc_name"])
    for row in rows:
        bgc_info[row["bgc_name"]] = {
            "organism": "TODO",
            "records": "TODO"
        }
    return bgc_info


def get_mibig_id_list(database):
    """returns a list of all bgc ids associated with MIBiG input files"""
    mibig_bgc_ids = [
        row["id"] for row in database.select(
            "bgc",
            "where dataset id = (",
            "select id from dataset",
            "where name = 'mibig'",
            ")",
            props=["id"])]
    return mibig_bgc_ids
