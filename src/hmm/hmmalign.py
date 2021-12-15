import os
import subprocess

from multiprocessing import Pool
from functools import partial

from src.utility.profiling import timeit
from src.pfam.stockholm import stockholm_parser


# @timeit
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
    