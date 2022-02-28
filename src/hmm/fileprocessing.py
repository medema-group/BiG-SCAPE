import logging
import os
import sys
from glob import glob

def get_cached_fasta_files(run):
    """Get all pre-cached fasta files"""
    file_path = os.path.join(run.directories.bgc_fasta, "*.fasta")
    fasta_files = glob(file_path)
    return set(fasta_files)

def get_cached_domain_fasta_files(run):
    """Get the cached domain fasta files (in output/cache/domains/*fasta)"""
    file_path = os.path.join(run.directories.domains, "*.fasta")
    fasta_files = glob(file_path)
    return set(fasta_files)

def get_searched_fasta_files(run, fasta_files):
    """From a set of fasta files, returns the list of fasta files that have an associated .pfd and
    .pfs file
    """
    # ignore this step and just return all of them if we are forcing hmmscan
    if run.options.force_hmmscan:
        logging.info(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
        return set()

    # assemble set
    already_done = set()
    for fasta in fasta_files:
        outputbase = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
        pfs_file_path = os.path.join(run.directories.pfs, outputbase + '.pfs')
        pfd_file_path = os.path.join(run.directories.pfd, outputbase + '.pfd')

        # ignore directories
        if os.path.isfile(pfs_file_path) and os.path.isfile(pfd_file_path):
            already_done.add(fasta)
    return already_done

def get_fasta_files_to_process(fasta_files, already_done):
    """Returns the subset of fasta_files which do not occur in already_done"""
    task_set = fasta_files - already_done
    if len(task_set) == 0:
        logging.info(" All fasta files had already been processed")
    elif len(already_done) > 0:
        if len(task_set) < 20:
            tasks = ", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in task_set)
            logging.warning(" The following NEW fasta file(s) will be processed: %s", tasks)
        else:
            logging.warning(" %d NEW fasta files will be processed", len(task_set))
    else:
        logging.info(" Predicting domains for %d fasta files", len(fasta_files))

    return list(task_set)


def check_fasta_files(run, cluster_base_names, all_fasta_files):
    """verify that all clusters have a corresponding fasta file in cache"""
    # fastaFiles: all the fasta files that should be there
    # (i.e. correspond to the input files)
    existing_fasta_files = set()
    for name in cluster_base_names:
        corresponding_name = os.path.join(run.directories.bgc_fasta, name+".fasta")
        existing_fasta_files.add(corresponding_name)

    # fastaBases: the actual fasta files we have that correspond to the input
    fasta_bases = all_fasta_files.intersection(existing_fasta_files)

    # Verify that all input files had their fasta sequences extracted
    if len(existing_fasta_files - fasta_bases) > 0:
        logging.error("The following files did NOT have their fasta sequences extracted: ")
        unextracted_files = existing_fasta_files - fasta_bases
        for unextracted_file in unextracted_files:
            logging.error(unextracted_file)
        sys.exit(1)


def get_cached_pfd_files(run):
    """Get the set of cached pfd files"""
    file_path = os.path.join(run.directories.pfd, "*.pfd")
    pfd_files = glob(file_path)
    return set(pfd_files)


def check_pfd_files(run, cluster_base_names):
    """Check the validity of the pfd files"""
    all_pfd_files = get_cached_pfd_files(run)

    # pfdFiles: all pfd files corresponding to the input files
    # (some input files could've been removed due to not having predicted domains)
    pfd_files = set()
    for name in cluster_base_names:
        pfd_files.add(os.path.join(run.directories.pfd, name+".pfd"))

    # pfdBases: the actual set of input files that have pfd files
    pfd_bases = all_pfd_files.intersection(pfd_files)

    # verify previous step.
    # All BGCs without predicted domains should no longer be in baseNames
    if len(pfd_files - pfd_bases) > 0:
        logging.error("The following files did NOT have their domtable files processed:")
        unprocessed_domtable_files = pfd_files - pfd_bases
        for unprocessed_domtable_file in unprocessed_domtable_files:
            logging.error(unprocessed_domtable_file)
        sys.exit(1)
