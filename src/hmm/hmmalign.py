import logging
import os
import subprocess
import sys

from multiprocessing import Pool
from functools import partial
from glob import glob

from src.utility import get_fasta_keys, fasta_parser
from src.pfam.stockholm import stockholm_parser


def do_multiple_align(run, try_resume):
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
        launch_hmmalign(run.options.cores, domain_sequence_list, run.directories.pfam,
                        run.options.verbose)

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


def launch_hmmalign(cores, domain_sequence_list, pfam_dir, verbose):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    pool = Pool(cores, maxtasksperchild=32)
    partial_func = partial(run_hmmalign, pfam_dir=pfam_dir, verbose=verbose)
    pool.map(partial_func, domain_sequence_list)
    pool.close()
    pool.join()


def run_hmmalign(domain_file, pfam_dir, verbose):
    #domain_file already contains the full path, with the file extension
    domain_base = domain_file.split(os.sep)[-1][:-6]
    hmmfetch_pars = ["hmmfetch", os.path.join(pfam_dir, "Pfam-A.hmm.h3m"), domain_base]
    proc_hmmfetch = subprocess.Popen(hmmfetch_pars, stdout=subprocess.PIPE, shell=False)

    domain_file_stk = domain_file[:-6]+".stk"
    hmmalign_pars = ["hmmalign", "-o", domain_file_stk, "-", domain_file]
    proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout,
                                     stdout=subprocess.PIPE, shell=False)

    proc_hmmfetch.stdout.close()
    proc_hmmalign.communicate()[0]
    proc_hmmfetch.wait()

    logging.debug(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))

    stockholm_parser(domain_file_stk)
    #SeqIO.convert(domain_file_stk, "stockholm", domain_file[:-6]+".algn", "fasta")


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
