from io import StringIO
import io
import logging
import math
import os

from multiprocessing import Queue, Process

import pyhmmer

from src.data import Database
from src.data.bgc import BGC
from src.data.hmm import from_accession
from src.data.hsp import get_hsp_id, insert_hsp, insert_hsp_alignment
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
                    accession = domain.alignment.hmm_accession.decode()
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


def run_pyhmmer_worker(input_queue, output_queue, profiles, pipeline, database: Database):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        bgc_id = input_queue.get(True)
        if bgc_id is None:
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

        output_queue.put((bgc_id, hsps))


def run_pyhmmer(run, database: Database, ids_todo):
    """Scan a list of fastas using pyhmmer scan

    inputs:
        run: run details for this execution of BiG-SCAPE
        task_set: a list of fasta file paths

    returns:
        a list of hits
    """
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    # get hmm profiles
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        profiles = list(hmm_file.optimized_profiles())


    pipeline = pyhmmer.plan7.Pipeline(pyhmmer.easel.Alphabet.amino(), Z=len(profiles), bit_cutoffs="trusted")

    num_processes = run.options.cores

    working_q = Queue(num_processes)

    num_tasks = len(ids_todo)

    output_q = Queue(num_tasks)

    processes = []
    for thread_num in range(num_processes):
        thread_name = f"distance_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)
        new_process = Process(target=run_pyhmmer_worker, args=(working_q, output_q, profiles, pipeline, database))
        processes.append(new_process)
        new_process.start()

    id_idx = 0
    ids_done = 0

    hsps = []

    while True:
        all_tasks_put = id_idx == num_tasks
        all_tasks_done = ids_done == num_tasks

        if all_tasks_put and all_tasks_done:
            break
        
        # if not all_tasks_put:
        if not working_q.full() and not all_tasks_put:
            bgc_id = ids_todo[id_idx]
            working_q.put(bgc_id)
            id_idx += 1
            if not working_q.full():
                continue
           
        if not output_q.empty():
            bgc_id, task_hsps = output_q.get()

            filtered_hsps = filter_overlap(task_hsps, run.options.domain_overlap_cutoff)

            for hsp in filtered_hsps:
                
                serial_nr = hsp[0]
                cds_id = hsp[1]
                hmm_id = hsp[2]
                bitscore = hsp[3]
                # insert hsp
                insert_hsp(database, serial_nr, cds_id, hmm_id, bitscore)
                hsps.append(hsp)

                # commit every 500 hsps
                if len(hsps) % 500 == 0:
                    database.commit_inserts()

            # update bgc status when done
            update_bgc_status(database, bgc_id, 2)

            ids_done += 1

            # commit every 500 bgcs also
            if ids_done % 500 == 0:
                database.commit_inserts()

            # print progress every 10%
            if ids_done % math.ceil(num_tasks / 10) == 0:
                percent_done = ids_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, ids_done, num_tasks)
        # logging.info("adding result (now %d)", len(network_matrix))

    # just making sure
    database.commit_inserts()

    # insert alignments. Has to be done after inserts because only then are ids available
    for idx, hsp in enumerate(hsps):
        serial_nr, cds_id, hmm_id, bitscore, env_start, env_end, model_start, model_end, cds_start, cds_end, model_gaps, cds_gaps = hsp

        # get hsp id
        hsp_id = get_hsp_id(database, serial_nr, cds_id, hmm_id)
        if hsp_id is None:
            logging.error("Could not find hsp_id associated with newly added hsp")


        # insert hsp_alignment
        insert_hsp_alignment(database, hsp_id, env_start, env_end, model_start, model_end, model_gaps, cds_start, cds_end, cds_gaps)

        # commit every 500 rows
        if idx % 500 == 0:
            database.commit_inserts()

    database.commit_inserts()

    # clean up threads
    for thread_num in range(num_processes * 2):
        working_q.put(None)

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)
