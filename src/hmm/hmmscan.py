import logging
import math
import os

from multiprocessing import Queue, Process

import pyhmmer


from src.pfam.misc import check_overlap


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
    cutoffs = dict()
    for profile in profiles:
        cutoffs[profile.accession.decode()] = profile.cutoffs.trusted

        try:
            search_result = pipeline.search_hmm(profile, sequences)
        except TypeError:
            logging.warning("    parsing %s threw TypeError. Ignoring...", base_name)
            return

        for hit in search_result:
            if hit.is_included():
                for domain in hit.domains:
                    accession = domain.alignment.hmm_accession.decode()
                    if domain.score > cutoffs[accession][1]:
                        yield domain

def run_pyhmmer_worker(input_queue, output_queue, profiles, pipeline):
    """worker for the run_pyhmmer method"""
    alphabet = pyhmmer.easel.Alphabet.amino()
    while True:
        base_name, fasta_file = input_queue.get(True)
        if base_name is None:
            break

        with pyhmmer.easel.SequenceFile(fasta_file) as seq_file:
            seq_file.set_digital(alphabet)
            try:
                sequences = list(seq_file)
            except ValueError:
                logging.warning("  Parsing %s threw ValueError.", base_name)
                logging.warning("  This means something in the source fasta file is wrong.")
                logging.warning("  A known case is a protein starting with a - instead of M.")
                logging.warning("  For the time being, this FASTA file will be skipped")
                output_queue.put((base_name, []))
                continue
            domains = pyhmmer_search_hmm(base_name, profiles, sequences, pipeline)
            del sequences

        pfd_matrix = []
        for domain in domains:
            score = f"{domain.score:.2f}"
            env_from = str(domain.env_from)
            env_to = str(domain.env_to)
            pfam_id = domain.alignment.hmm_accession.decode()
            domain_name = domain.alignment.hmm_name.decode()

            header = domain.alignment.target_name.decode()
            header_list = header.split(":")
            try:
                gene_id = header_list[header_list.index("gid")+1]
            except ValueError:
                logging.warning("  No gene ID in %s", base_name)
                gene_id = ''
            gene_start = header_list[header_list.index("loc")+1]
            gene_end = header_list[header_list.index("loc")+2]

            pfd_row = [base_name]
            pfd_row.append(score)
            pfd_row.append(gene_id)
            pfd_row.append(env_from)
            pfd_row.append(env_to)
            pfd_row.append(pfam_id)
            pfd_row.append(domain_name)
            pfd_row.append(gene_start)
            pfd_row.append(gene_end)
            pfd_row.append(header)

            pfd_matrix.append(pfd_row)
            del pfd_row

        output_queue.put((base_name, pfd_matrix))


def run_pyhmmer(run, fasta_files_to_process, genbank_dict, clusters, base_names, mibig_set):
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
    # num_processes = 8

    working_q = Queue(num_processes)

    num_tasks = len(fasta_files_to_process)

    output_q = Queue(num_tasks)

    processes = []
    for thread_num in range(num_processes):
        thread_name = f"distance_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)
        new_process = Process(target=run_pyhmmer_worker, args=(working_q, output_q, profiles, pipeline))
        processes.append(new_process)
        new_process.start()

    pfd_matrix = []
    fasta_idx = 0
    fastas_done = 0

    while True:
        all_tasks_put = fasta_idx == num_tasks
        all_tasks_done = fastas_done == num_tasks

        if all_tasks_put and all_tasks_done:
            break

        if not working_q.full() and not all_tasks_put:
            fasta_file = fasta_files_to_process[fasta_idx]
            base_name = ".".join(fasta_file.split(os.sep)[-1].split(".")[:-1])
            working_q.put((base_name, fasta_file))
            fasta_idx += 1
            if not working_q.full():
                continue

        if not output_q.empty():
            base_name, pfd_matrix = output_q.get()

            if len(pfd_matrix) == 0:
                # there aren't any domains in this BGC
                # delete from all data structures
                logging.info("  No domains were found in %s. Removing it from further analysis", base_name)

                pfdoutput = os.path.join(run.directories.pfd, base_name + ".pfd")
                clusters.remove(base_name)
                base_names.remove(base_name)
                del genbank_dict[base_name]
                if base_name in mibig_set:
                    mibig_set.remove(base_name)
            else:

                # check_overlap also sorts the filtered_matrix results and removes
                # overlapping domains, keeping the highest scoring one
                filtered_matrix, domains = check_overlap(pfd_matrix, run.options.domain_overlap_cutoff)

                # Save list of domains per BGC
                pfsoutput = os.path.join(run.directories.pfs, base_name + ".pfs")
                with open(pfsoutput, 'w', encoding="UTF-8") as pfs_handle:
                    pfs_handle.write(" ".join(domains))

                # Save more complete information of each domain per BGC
                pfdoutput = os.path.join(run.directories.pfd, base_name + ".pfd")
                with open(pfdoutput, 'w', encoding="UTF-8") as pfd_handle:
                    write_pfd(pfd_handle, filtered_matrix)

            fastas_done += 1

            # print progress every 10%
            if fastas_done % math.ceil(num_tasks / 10) == 0:
                percent_done = fastas_done / num_tasks * 100
                logging.info("  %d%% (%d/%d)", percent_done, fastas_done, num_tasks)
            # logging.info("adding result (now %d)", len(network_matrix))

    # clean up threads
    for thread_num in range(num_processes):
        working_q.put((None, None))

    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)
