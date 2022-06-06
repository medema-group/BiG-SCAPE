import logging
import math
import multiprocessing
import os

from multiprocessing.connection import Connection, wait
from multiprocessing import Pipe, Process

import pyhmmer

from src.data.functions import get_bgc_name_by_id, remove_bgc
from src.data import Database
from src.data.bgc import BGC
from src.data.hmm import from_accession
from src.data.hsp import get_hsp_id, insert_feature, insert_hsp, insert_hsp_alignment
from src.data.status import update_bgc_status


def write_pfd(pfd_handle, matrix):
    """Writes all rows in a pfd matrix to the file given by pfd_handle"""
    for row in matrix:
        row = "\t".join(row)
        pfd_handle.write(row+"\n")

    pfd_handle.close()

def pyhmmer_hmmpress(options):
    """Presses the Pfam-A.hmm file into optimized files for further analysis

    Inputs:
        run: run details for this execution of BiG-SCAPE
    """
    hmm_file_path = os.path.join(options.pfam_dir, "Pfam-A.hmm")
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        if hmm_file.is_pressed():
            logging.info(" PFAM Hmm file was already pressed.")
            return

        logging.info(" Pressing Pfam HMM file...")
        pyhmmer.hmmer.hmmpress(hmm_file, os.path.join(options.pfam_dir, "Pfam-A.hmm"))

def pyhmmer_search_hmm(base_name, profiles, sequences, pipeline):
    """scan an iterable of pyhmmer sequences using an optimized profiles

    Inputs:
        profiles: optimized profiles from HMMFile.optimized_profiles()
        pipeline: pipeline object for HMM searching
        sequences: iterable of sequences

    Yields:
        A list of domains, similar to what hmmsearch would return
    """
    for profile in profiles:
        try:
            search_result = pipeline.search_hmm(profile, sequences)
        except TypeError:
            logging.warning("    parsing %s threw TypeError. Ignoring...", base_name)
            return

        for hit in search_result:
            if hit.is_included():
                for domain in hit.domains:
                    if domain.score > profile.cutoffs.trusted[1]:
                        yield domain


def get_cds_gaps(cds_sequence):
    "gets a comma-separated list of gap locations in the cds sequence"
    gaps = []
    for idx, aa in enumerate(cds_sequence):
        if aa == "-":
            gaps.append(idx)
    return ",".join(map(str, gaps))

def get_hmm_gaps(hmm_sequence):
    "gets a comma-separated list of gap locations in the model sequence"
    gaps = []
    for idx, aa in enumerate(hmm_sequence):
        if aa == ".":
            gaps.append(idx)
    return ",".join(map(str, gaps))

def no_overlap(a_start, a_end, b_start, b_end):
    """Return True if there is no overlap between two regions"""
    if a_start < b_start and a_end < b_start:
        return True
    elif a_start > b_end and a_end > b_end:
        return True
    else:
        return False

def len_overlap(a_start, a_end, b_start, b_end):
    """Returns the length of an overlapping sequence"""

    if a_start < b_start:
        cor1 = a_start
    else:
        cor1 = b_start

    if a_end > b_end:
        cor2 = a_end
    else:
        cor2 = b_end

    total_region = cor2 - cor1
    sum_len = (a_end - a_start) + (b_end - b_start)

    return sum_len - total_region

def filter_overlap(hsps, overlap_cutoff):
    """Check if domains overlap for a certain overlap_cutoff.
     If so, remove the domain(s) with the lower score."""

    delete_list = []
    for i in range(len(hsps)-1):
        for j in range(i+1, len(hsps)):
            row1 = hsps[i]
            row2 = hsps[j]

            # using env coords
            a_start = row1[4]
            a_end = row1[5]
            b_start = row2[4]
            b_end = row2[5]

            #check if we are the same CDS
            if row1[1] == row2[1]:
                #check if there is overlap between the domains
                if not no_overlap(a_start, a_end, b_start, b_end):
                    overlapping_aminoacids = len_overlap(a_start, a_end, b_start, b_end)
                    overlap_perc_loc1 = overlapping_aminoacids / (a_end - a_start)
                    overlap_perc_loc2 = overlapping_aminoacids / (b_end - b_start)
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[3]) >= float(row2[3]): #see which has a better score
                            delete_list.append(row2)
                        elif float(row1[3]) < float(row2[3]):
                            delete_list.append(row1)

    for lst in delete_list:
        try:
            hsps.remove(lst)
        except ValueError:
            pass
    return hsps

def rank_normalize_hsps(hsps, top_k):
    """Rank normalizes a set of hsps per cds
    """
    sorted_hsps = sorted(hsps, key = lambda elem: elem[3], reverse=True)
    result_hsps = []
    cds_hsp_count = dict()
    for hsp in sorted_hsps:
        serial_nr = hsp[0]
        cds_id = hsp[1]
        hmm_id = hsp[2]
        value = int(hsp[3])
        if cds_id not in cds_hsp_count:
            cds_hsp_count[cds_id] = 0
        if top_k > 0 and top_k < cds_hsp_count[cds_id] + 1:
            continue
        value = int(255 - int((255 / top_k) * cds_hsp_count[cds_id]))
        cds_hsp_count[cds_id] += 1
        result_hsps.append([serial_nr, cds_id, hmm_id, value])
    return result_hsps


def run_pyhmmer_worker(connection: Connection, profiles, pipeline, database: Database):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        bgc_id = connection.recv()
        if bgc_id is None:
            connection.close()
            break

        base_name = BGC.get_bgc_base_name(bgc_id, database)
        bgc_cds_list = BGC.get_all_cds([bgc_id], database)

        sequences = []
        for cds_row in bgc_cds_list:
            accession = BGC.CDS.gen_accession(base_name, cds_row)
            ds = pyhmmer.easel.TextSequence(accession=accession.encode(), name=str(cds_row["id"]).encode(), sequence=cds_row["aa_seq"]).digitize(alphabet)
            sequences.append(ds)

        domains = pyhmmer_search_hmm(accession, profiles, sequences, pipeline)

        hsps = list()
        for idx, domain in enumerate(domains):
            domain: pyhmmer.plan7.Domain
            # cds id
            cds_id = int(domain.alignment.target_name.decode())

            # hmm id
            hmm_accession = domain.alignment.hmm_accession.decode()
            # only happens on subpfams
            if hmm_accession == "" or hmm_accession is None:
                hmm_accession = domain.alignment.hmm_name.decode()
            hmm = from_accession(database, hmm_accession)
            hmm_id = hmm["id"]

            # score
            bitscore = domain.score

            # env coords
            env_start = domain.env_from
            env_end = domain.env_to

            # model coords
            model_start = domain.alignment.hmm_from
            model_end = domain.alignment.hmm_to

            # cds coords
            cds_start = domain.alignment.target_from
            cds_end = domain.alignment.target_to

            # gaps
            model_gaps = get_hmm_gaps(domain.alignment.hmm_sequence)
            cds_gaps = get_cds_gaps(domain.alignment.target_sequence)

            hsps.append((idx, cds_id, hmm_id, bitscore, env_start, env_end, model_start, model_end, cds_start, cds_end, model_gaps, cds_gaps))

        connection.send((bgc_id, hsps))

def run_pyhmmer_pfam(run, database: Database, ids_todo):
    """Runs the pyhmmer pipeline using the pfam hmm"""
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    # get hmm profiles
    # pfam
    profiles = list()
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        profiles.extend(list(hmm_file.optimized_profiles()))
        logging.info("Found %d hmm profiles", len(profiles))

    run_pyhmmer(run, database, "hsp", ids_todo, profiles)


def run_pyhmmer(
    run,
    database: Database,
    hsp_table,
    ids_todo,
    profiles,
    use_filter_overlap=True,
    insert_alignments=True,
    generate_features=False
):
    """Scan a list of fastas using pyhmmer scan

    inputs:
        run: run details for this execution of BiG-SCAPE
        task_set: a list of fasta file paths

    returns:
        a list of hits
    """
    pipeline = pyhmmer.plan7.Pipeline(pyhmmer.easel.Alphabet.amino(), Z=len(profiles), bit_cutoffs="trusted")

    num_processes = run.options.cores

    connections = list()

    num_tasks = len(ids_todo)

    id_idx = 0
    ids_done = 0

    hsps = []
    # keep track of bgs with no domains
    bgc_no_domains = []

    processes = []
    for thread_num in range(num_processes):
        # if there are fewer tasks than there are threads, don't start a new thread
        # this is only a concern for < cores * 4 tasks
        if id_idx >= num_tasks:
            continue
        # generate thread name
        thread_name = f"hmmscan_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)
        
        # generate connection
        main_connection, worker_connection = Pipe(True)
        connections.append(main_connection)

        # create and start new process
        new_process = Process(target=run_pyhmmer_worker, args=(worker_connection, profiles, pipeline, database))
        processes.append(new_process)
        new_process.start()

        # send starting data
        bgc_id = ids_todo[id_idx]
        main_connection.send(bgc_id)
        id_idx += 1

    while True:
        all_tasks_put = id_idx == num_tasks
        all_tasks_done = ids_done == num_tasks

        # break if everything is done
        if all_tasks_put and all_tasks_done:
            break

        available_connections = wait(connections)
        for connection in available_connections:
            connection: Connection

            # receive object
            bgc_id, task_hsps = connection.recv()

            # add new tasks first
            if all_tasks_put:
                # if done close this process
                connection.send(None)
                connection.close()
                connections.remove(connection)
            else:
                # else add a new task
                next_bgc_id = ids_todo[id_idx]
                connection.send(next_bgc_id)
                id_idx += 1
                all_tasks_put = id_idx == num_tasks

            # process results
            result_hsps: list = task_hsps

            if use_filter_overlap:
                result_hsps = filter_overlap(result_hsps, run.options.domain_overlap_cutoff)

            if len(result_hsps) == 0:
                bgc_no_domains.append(bgc_id)

            for idx, hsp in enumerate(result_hsps):
                serial_nr = hsp[0]
                cds_id = hsp[1]
                hmm_id = hsp[2]
                bitscore = hsp[3]
                # insert hsp
                insert_hsp(database, hsp_table, serial_nr, cds_id, hmm_id, bitscore)
                hsps.append(hsp)

                # commit every 500 hsps
                if len(hsps) % 500 == 0:
                    database.commit_inserts()

            if generate_features:
                # order by bitscore if rank_normalize is true
                rank_normalized_hsps = rank_normalize_hsps(result_hsps, 3)

                for idx, hsp in enumerate(rank_normalized_hsps):
                    # insert hsp
                    insert_feature(database, hsp)

                    # commit every 500 hsps
                    if len(hsps) % 500 == 0:
                        database.commit_inserts()


            # update bgc status when done
            update_bgc_status(database, bgc_id, 2)

            ids_done += 1
            all_tasks_done = ids_done == num_tasks

            # commit every 500 bgcs also
            if ids_done % 500 == 0:
                database.commit_inserts()

            # print progress every 10%
            if ids_done % math.ceil(num_tasks / 10) == 0:
                percent_done = ids_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, ids_done, num_tasks)

    # just making sure
    database.commit_inserts()

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)

    # clean up bgcs that didn't have any domains
    for bgc_id in bgc_no_domains:
        bgc_name = get_bgc_name_by_id(database, bgc_id)
        logging.warning("    BGC %s has no domains and will be removed from \
            further analysis", bgc_name)
        remove_bgc(database, bgc_id)

    # insert alignments. Has to be done after inserts because only then are ids available
    if not insert_alignments:
        return

    for idx, hsp in enumerate(hsps):
        serial_nr, cds_id, hmm_id, bitscore, env_start, env_end, model_start, model_end, cds_start, cds_end, model_gaps, cds_gaps = hsp

        # get hsp id
        hsp_id = get_hsp_id(database, serial_nr, cds_id, hmm_id)
        if hsp_id is None:
            logging.error("    Could not find hsp_id associated with newly added hsp")

        # insert hsp_alignment
        insert_hsp_alignment(database, hsp_id, env_start, env_end, model_start, model_end, model_gaps, cds_start, cds_end, cds_gaps)

        # commit every 500 rows
        if idx % 500 == 0:
            database.commit_inserts()

    database.commit_inserts()
