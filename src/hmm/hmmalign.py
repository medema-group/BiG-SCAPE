import logging
import math
import os
import subprocess
import sys
import pyhmmer

from multiprocessing import Pool, Process, Queue
from functools import partial
from glob import glob

from src.utility import get_fasta_keys, fasta_parser
from src.pfam.stockholm import stockholm_parser


def do_multiple_align(run, try_resume):
    """Perform multiple alignment of domain sequences

    inputs:
        run: run details for this execution of BiG-SCAPE
        try_resume: boolean indicating if BiG-SCAPE detected an interrupted alignment
    """
    logging.info("Performing multiple alignment of domain sequences")
    # obtain all fasta files with domain sequences
    domain_sequence_list = set(glob(os.path.join(run.directories.domains, "*.fasta")))

    # compare with .algn set of files. Maybe resuming is possible if
    # no new sequences were added
    if try_resume:
        temp_aligned = set(glob(os.path.join(run.directories.domains, "*.algn")))

        if len(temp_aligned) > 0:
            logging.info(" Found domain fasta files without corresponding alignments")

            for algn_file in temp_aligned:
                if os.path.getsize(algn_file) > 0:
                    domain_sequence_list.remove(algn_file[:-5]+".fasta")

        temp_aligned.clear()

    # Try to further reduce the set of domain fastas that need alignment
    sequence_tag_list = set()
    header_list = []
    domain_sequence_list_temp = domain_sequence_list.copy()
    for domain_file in domain_sequence_list_temp:
        domain_name = ".".join(domain_file.split(os.sep)[-1].split(".")[:-1])

        # fill fasta_dict...
        with open(domain_file, "r") as fasta_handle:
            header_list = get_fasta_keys(fasta_handle)

        # Get the BGC name from the sequence tag. The form of the tag is:
        # >BGCXXXXXXX_BGCXXXXXXX_ORF25:gid...
        sequence_tag_list = set(s.split("_ORF")[0] for s in header_list)

        # ...to find out how many sequences do we actually have
        if len(sequence_tag_list) == 1:
            # avoid multiple alignment if the domains all belong to the same BGC
            domain_sequence_list.remove(domain_file)

            logging.debug(" Skipping Multiple Alignment for %s \
                   (appears only in one BGC)", domain_name)

    sequence_tag_list.clear()
    del header_list[:]

    domain_sequence_list_temp.clear()

    # Do the multiple alignment
    stop_flag = False
    if len(domain_sequence_list) > 0:
        logging.info("\n Using hmmalign")
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

        launch_hmmalign(run, domain_sequence_list, profile_dictionary)

        # verify all tasks were completed by checking existance of alignment files
        for domain_file in domain_sequence_list:
            if not os.path.isfile(domain_file[:-6]+".algn"):
                logging.error("   %s.algn could not be found \
                    (possible issue with aligner).", domain_file[:-6])
                stop_flag = True
        if stop_flag:
            sys.exit(1)
    else:
        logging.info(" No domain fasta files found to align")

def launch_hmmalign_worker(input_queue, output_queue, profile_dict):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        domain_fasta, temp = input_queue.get(True)
        if domain_fasta is None:
            break

        base_name = ".".join(domain_fasta.split(os.sep)[-1].split(".")[:-1])
        profile = profile_dict[base_name]
        with pyhmmer.easel.SequenceFile(domain_fasta) as seq_file:
            seq_file.set_digital(alphabet)
            sequences = list(seq_file)

            # do the alignment
            msa = pyhmmer.hmmer.hmmalign(profile, sequences)

            # write stk TODO: remove, write stuff directly
            domain_file_stk = domain_fasta[:-6]+".stk"
            with open(domain_file_stk, "wb") as stk_file:
                msa.write(stk_file, "stockholm")

            # parse stockholm file
            stockholm_parser(domain_file_stk)

            # done, write something to output
            output_queue.put((domain_fasta, None))


def launch_hmmalign(run, domain_sequence_set, profile_dict):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    # prepare data
    num_processes = run.options.cores
    # num_processes = 8

    working_q = Queue(num_processes)

    domain_sequence_list = list(domain_sequence_set)
    num_tasks = len(domain_sequence_list)

    output_q = Queue(num_tasks)

    processes = []
    for thread_num in range(num_processes):
        thread_name = f"distance_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)
        new_process = Process(target=launch_hmmalign_worker, args=(working_q, output_q, profile_dict))
        processes.append(new_process)
        new_process.start()

    fasta_idx = 0
    fastas_done = 0

    while True:
        all_tasks_put = fasta_idx == num_tasks
        all_tasks_done = fastas_done == num_tasks

        if all_tasks_put and all_tasks_done:
            break

        if not working_q.full() and not all_tasks_put:
            fasta_file = domain_sequence_list[fasta_idx]
            working_q.put((fasta_file, None))
            fasta_idx += 1
            if not working_q.full():
                continue

        if not output_q.empty():
            fasta_file, temp = output_q.get()
        
            fastas_done += 1

            # print progress every 10%
            if fastas_done % math.ceil(num_tasks / 10) == 0:
                percent_done = fastas_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, fastas_done, num_tasks)

    # clean up threads
    for thread_num in range(num_processes):
        working_q.put((None, None))

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)


def read_aligned_files(run):
    aligned_domain_seqs = {} # Key: specific domain sequence label. Item: aligned sequence
    aligned_files_list = glob(os.path.join(run.directories.domains, "*.algn"))
    if len(aligned_files_list) == 0:
        logging.error("No aligned sequences found in the domain folder (run without the --skip_ma \
                 parameter or point to the correct output folder)")
        sys.exit(1)
    for aligned_file in aligned_files_list:
        with open(aligned_file, "r") as aligned_file_handle:
            fasta_dict = fasta_parser(aligned_file_handle)
            for header in fasta_dict:
                aligned_domain_seqs[header] = fasta_dict[header]
    return aligned_domain_seqs
