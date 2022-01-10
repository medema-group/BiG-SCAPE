import os
import subprocess
import sys

from multiprocessing import Pool
from functools import partial
from glob import glob

from src.utility import get_fasta_keys, fasta_parser
from src.pfam.stockholm import stockholm_parser


def do_multiple_align(run, try_resume):
    print("Performing multiple alignment of domain sequences")
    # obtain all fasta files with domain sequences
    DOMAIN_SEQUENCE_LIST = set(glob(os.path.join(run.directories.domains, "*.fasta")))

    # compare with .algn set of files. Maybe resuming is possible if
    # no new sequences were added
    if try_resume:
        TEMP_ALIGNED = set(glob(os.path.join(run.directories.domains, "*.algn")))

        if len(TEMP_ALIGNED) > 0:
            print(" Found domain fasta files without corresponding alignments")

            for a in TEMP_ALIGNED:
                if os.path.getsize(a) > 0:
                    DOMAIN_SEQUENCE_LIST.remove(a[:-5]+".fasta")

        TEMP_ALIGNED.clear()

    # Try to further reduce the set of domain fastas that need alignment
    SEQUENCE_TAG_LIST = set()
    HEADER_LIST = []
    DOMAIN_SEQUENCE_LIST_TEMP = DOMAIN_SEQUENCE_LIST.copy()
    for domain_file in DOMAIN_SEQUENCE_LIST_TEMP:
        domain_name = ".".join(domain_file.split(os.sep)[-1].split(".")[:-1])

        # fill fasta_dict...
        with open(domain_file, "r") as fasta_handle:
            HEADER_LIST = get_fasta_keys(fasta_handle)

        # Get the BGC name from the sequence tag. The form of the tag is:
        # >BGCXXXXXXX_BGCXXXXXXX_ORF25:gid...
        SEQUENCE_TAG_LIST = set(s.split("_ORF")[0] for s in HEADER_LIST)

        # ...to find out how many sequences do we actually have
        if len(SEQUENCE_TAG_LIST) == 1:
            # avoid multiple alignment if the domains all belong to the same BGC
            DOMAIN_SEQUENCE_LIST.remove(domain_file)
            if run.options.verbose:
                print(" Skipping Multiple Alignment for {} \
                    (appears only in one BGC)".format(domain_name))

    SEQUENCE_TAG_LIST.clear()
    del HEADER_LIST[:]

    DOMAIN_SEQUENCE_LIST_TEMP.clear()

    # Do the multiple alignment
    STOP_FLAG = False
    if len(DOMAIN_SEQUENCE_LIST) > 0:
        print("\n Using hmmalign")
        launch_hmmalign(run.options.cores, DOMAIN_SEQUENCE_LIST, run.directories.pfam,
                            run.options.verbose)

        # verify all tasks were completed by checking existance of alignment files
        for domain_file in DOMAIN_SEQUENCE_LIST:
            if not os.path.isfile(domain_file[:-6]+".algn"):
                print("   ERROR, {}.algn could not be found \
                    (possible issue with aligner).".format(domain_file[:-6]))
                STOP_FLAG = True
        if STOP_FLAG:
            sys.exit()
    else:
        print(" No domain fasta files found to align")


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
    hmmfetch_pars = ["hmmfetch", os.path.join(pfam_dir,"Pfam-A.hmm.h3m"), domain_base]
    proc_hmmfetch = subprocess.Popen(hmmfetch_pars, stdout=subprocess.PIPE, shell=False)
    
    domain_file_stk = domain_file[:-6]+".stk"
    hmmalign_pars = ["hmmalign", "-o", domain_file_stk, "-", domain_file]
    proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout, stdout=subprocess.PIPE, shell=False)
    
    proc_hmmfetch.stdout.close()
    proc_hmmalign.communicate()[0]
    proc_hmmfetch.wait()
    
    if verbose:
        print(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))
    
    stockholm_parser(domain_file_stk)
    #SeqIO.convert(domain_file_stk, "stockholm", domain_file[:-6]+".algn", "fasta")


def read_aligned_files(run):
    aligned_domain_seqs = {} # Key: specific domain sequence label. Item: aligned sequence
    ALIGNED_FILES_LIST = glob(os.path.join(run.directories.domains, "*.algn"))
    if len(ALIGNED_FILES_LIST) == 0:
        sys.exit("No aligned sequences found in the domain folder (run without the --skip_ma parameter or point to the correct output folder)")
    for aligned_file in ALIGNED_FILES_LIST:
        with open(aligned_file, "r") as aligned_file_handle:
            fasta_dict = fasta_parser(aligned_file_handle)
            for header in fasta_dict:
                aligned_domain_seqs[header] = fasta_dict[header]
    return aligned_domain_seqs
