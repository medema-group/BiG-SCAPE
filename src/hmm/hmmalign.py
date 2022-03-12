import logging
import math
import os
import string
import subprocess
import sys
import pyhmmer

from multiprocessing import Pool, Process, Queue
from functools import partial
from glob import glob
from src.data.cds import gen_accession, get_cds_rows
from src.data.database import Database
from src.data.hsp import get_multiple_align_hsps
from src.data.msa import insert_msa
from src.data.status import update_bgc_status

from src.utility import get_fasta_keys, fasta_parser
from src.pfam.stockholm import stockholm_parser

def generate_task_list(hsp_rows):
    """Generates a task list suitable for launch_hmmalign from a list of rows from the hsp table
    This is in the form of a list with the following structure:
    [(accession, [data row])]
    """
    # first convert relevant information into a dictionary
    data_dict = dict()
    for row in hsp_rows:
        if row["accession"] not in data_dict:
            data_dict[row["accession"]] = []
        data_dict[row["accession"]].append(row)
    task_list = []
    for accession, rows in data_dict.items():
        task_list.append((accession, rows))

    return task_list

def do_multiple_align(run, database: Database, hsp_ids: list()):
    """Perform multiple alignment of high scoring protein domain sequences

    inputs:
        run: run details for this execution of BiG-SCAPE
        database: database object for this data set
        hsp_ids: List of high scoring protein domains to analyze
    """
    logging.info(" Performing multiple alignment of domain sequences")

    # remove any hmms with only one bgc
    hsp_rows = get_multiple_align_hsps(database)
    algn_task_list = generate_task_list(hsp_rows)

    # quitting early if there are no alignments left.
    if len(algn_task_list) == 0:
        logging.info(" No alignments actually need to be done")
        return


    stop_flag = False
    logging.info(" Using hmmalign")
    # load the hmm profiles
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    # get hmm profiles

    # load the hmm file
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        profiles = list(hmm_file)

    # we will pass the relevant profiles to the subprocesses later on
    # so we want to assemble a dictionary of profiles
    profile_dictionary = dict()
    for profile in profiles:
        profile_accession = profile.accession.decode()
        profile_dictionary[profile_accession] = profile
    # clear this data
    del profiles

    # hand off to the thread handler
    launch_hmmalign(run, algn_task_list, profile_dictionary, database)

def process_algn_string(algn_string: str):
    """removes any model gaps from the alignment string and returns a string with only cds gaps"""
    return algn_string.translate(str.maketrans('', '', string.ascii_lowercase + "."))

def launch_hmmalign_worker(input_queue, output_queue, profile_dict, database):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        accession, hsp_rows = input_queue.get(True)
        if accession is None:
            break

        # get profile
        profile = profile_dict[accession]

        # get hmm_id
        hmm_id = hsp_rows[0]["hmm_id"]


        # get cds sequences
        cds_ids = [hsp_row["cds_id"] for hsp_row in hsp_rows]
        cds_rows = get_cds_rows(database, cds_ids)

        sequences = list()
        for cds_row in cds_rows:
            ds = pyhmmer.easel.TextSequence(name=str(cds_row["id"]).encode(), sequence=cds_row["aa_seq"]).digitize(pyhmmer.easel.Alphabet.amino())
            sequences.append(ds)


        # do the alignment
        msa = pyhmmer.hmmer.hmmalign(profile, sequences)

        alignments = []
        for idx, alignment in enumerate(msa.alignment):
            cds_id = msa.sequences[idx].name.decode()
            algn_string = process_algn_string(alignment)
            alignments.append((cds_id, hmm_id, algn_string))

        # done, write something to output
        output_queue.put(alignments)


def launch_hmmalign(run, algn_task_list, profile_dict, database: Database):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    # prepare data
    num_processes = run.options.cores
    # num_processes = 8

    working_q = Queue(num_processes)

    num_tasks = len(algn_task_list)

    output_q = Queue(num_tasks)

    processes = []
    for thread_num in range(num_processes):
        thread_name = f"MSA_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)
        new_process = Process(target=launch_hmmalign_worker, args=(working_q, output_q, profile_dict, database))
        processes.append(new_process)
        new_process.start()

    task_idx = 0
    tasks_done = 0

    while True:
        all_tasks_put = task_idx == num_tasks
        all_tasks_done = tasks_done == num_tasks

        if all_tasks_put and all_tasks_done:
            break

        if not working_q.full() and not all_tasks_put:
            task = algn_task_list[task_idx]
            working_q.put(task)
            task_idx += 1
            if not working_q.full():
                continue

        if not output_q.empty():
            alignments = output_q.get()
            for alignment in alignments:
                cds_id, hmm_id, algn_string = alignment

                insert_msa(database, cds_id, hmm_id, algn_string)

            tasks_done += 1

            # commit every 500 results
            if tasks_done % 500 == 0:
                database.commit_inserts()

            # print progress every 10%
            if tasks_done % math.ceil(num_tasks / 10) == 0:
                percent_done = tasks_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, tasks_done, num_tasks)

    # commit changes to database
    database.commit_inserts()

    # clean up threads
    for thread_num in range(num_processes * 2):
        working_q.put((None, None))

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)

