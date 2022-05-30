import logging
import math
from multiprocessing.connection import Connection, wait
import os
import string
import pyhmmer

from multiprocessing import Pipe, Process
from src.data.database import Database
from src.data.hsp import get_multiple_align_hsps, get_hsp_cds
from src.data.msa import insert_msa

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

def launch_hmmalign_worker(connection: Connection, profile_dict, database):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        accession, hsp_rows = connection.recv()
        if accession is None:
            connection.close()
            break

        # get profile
        profile = profile_dict[accession]

        # get hmm_id
        hmm_id = hsp_rows[0]["hmm_id"]


        # get cds sequences
        cds_ids = [hsp_row["cds_id"] for hsp_row in hsp_rows]
        cds_rows = get_hsp_cds(database, cds_ids, hmm_id)

        sequences = list()
        for cds_row in cds_rows:
            cds_id = cds_row["cds_id"]
            env_start = cds_row["env_start"]
            env_end = cds_row["env_end"]
            name = f"{cds_id}_{env_start}_{env_end}".encode()
            sequence = sequence=cds_row["sequence"]
            ds = pyhmmer.easel.TextSequence(
                name = name,
                sequence = sequence
            ).digitize(pyhmmer.easel.Alphabet.amino())
            sequences.append(ds)


        # do the alignment
        msa = pyhmmer.hmmer.hmmalign(profile, sequences)

        alignments = []
        for idx, alignment in enumerate(msa.alignment):
            name = msa.sequences[idx].name.decode()
            cds_id, env_start, env_end = name.split("_")
            algn_string = process_algn_string(alignment)
            alignments.append((cds_id, hmm_id, env_start, env_end, algn_string))

        # done, write something to output
        connection.send(alignments)


def launch_hmmalign(run, algn_task_list, profile_dict, database: Database):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    # prepare data
    num_processes = run.options.cores

    connections = list()

    num_tasks = len(algn_task_list)

    task_idx = 0
    tasks_done = 0

    processes = []
    for thread_num in range(num_processes):
        # if there are fewer tasks than there are threads, don't start a new thread
        # this is only a concern for < cores * 4 tasks
        if task_idx >= num_tasks:
            continue
        # generate thread name
        thread_name = f"MSA_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)

        # generate connection
        main_connection, worker_connection = Pipe(True)
        connections.append(main_connection)

        # create and start new process
        new_process = Process(target=launch_hmmalign_worker, args=(worker_connection, profile_dict, database))
        processes.append(new_process)
        new_process.start()

        # send starting data
        task = algn_task_list[task_idx]
        main_connection.send(task)
        task_idx += 1

    while True:
        all_tasks_put = task_idx == num_tasks
        all_tasks_done = tasks_done == num_tasks

        if all_tasks_put and all_tasks_done:
            break

        available_connections = wait(connections)
        for connection in available_connections:
            connection: Connection

            # receive data
            alignments = connection.recv()

            # add new tasks first
            if all_tasks_put:
                # if done close this process
                connection.send((None, None))
                connection.close()
                connections.remove(connection)
            else:
                task = algn_task_list[task_idx]
                connection.send(task)
                task_idx += 1
                all_tasks_put = task_idx == num_tasks

            # process results
            for alignment in alignments:
                cds_id, hmm_id, env_start, env_end, algn_string = alignment

                insert_msa(database, cds_id, hmm_id, env_start, env_end, algn_string)

            tasks_done += 1
            all_tasks_done = tasks_done == num_tasks

            # commit every 500 results
            if tasks_done % 500 == 0:
                database.commit_inserts()

            # print progress every 10%
            if tasks_done % math.ceil(num_tasks / 10) == 0:
                percent_done = tasks_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, tasks_done, num_tasks)


    # commit changes to database
    database.commit_inserts()

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)

