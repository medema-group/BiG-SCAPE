import logging
import os
import sys
from glob import glob

def get_cached_fasta_files(run):
    # All available fasta files (could be more than it should if reusing output folder)
    file_path = os.path.join(run.directories.bgc_fasta, "*.fasta")
    fasta_files = glob(file_path)
    return set(fasta_files)

def get_cached_domain_fasta_files(run):
    """Get the cached domain fasta files (in output/cache/domains/*fasta)"""
    file_path = os.path.join(run.directories.domains, "*.fasta")
    fasta_files = glob(file_path)
    return set(fasta_files)


def check_fasta_files(run, cluster_base_names, all_fasta_files):
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

def get_fasta_files_to_process(run, fasta_files):
    """Gets the fasta files that do not have an associated pfd or pfs file"""
    if run.options.force_hmmscan:
        logging.info(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
        return fasta_files

    already_done = set()
    for fasta in fasta_files:
        outputbase = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
        pfs_file_path = os.path.join(run.directories.pfd, outputbase + '.pfs')
        pfd_file_path = os.path.join(run.directories.pfd, outputbase + '.pfd')

        # ignore directories
        if os.path.isfile(pfs_file_path) and os.path.isfile(pfd_file_path):
            already_done.add(fasta)


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

def get_cached_domtable_files(run):
    file_path = os.path.join(run.directories.domtable, "*.domtable")
    domtable_files = glob(file_path)
    return set(domtable_files)

def get_domtable_files_to_process(run, domtable_files):
    # return early if force_hmmscan is true
    if run.options.force_hmmscan:
        return domtable_files

    already_done = get_searched_fasta_files(run, domtable_files)

    # the domtable files to be processed is the difference between all, and those which were
    # already done
    unprocessed_domtable_files = domtable_files - already_done
    first_run = len(domtable_files) == 0

    log_unprocessed_domtable_files(unprocessed_domtable_files, first_run)

    return unprocessed_domtable_files

def get_searched_fasta_files(run, fasta_files):
    # otherwise find what still needs to be done
    already_done = set()
    for domtable in fasta_files:
        outputbase = ".".join(domtable.split(os.sep)[-1].split(".")[:-1])
        outputfile = os.path.join(run.directories.pfd, outputbase + '.pfd')
        if os.path.isfile(outputfile):
            already_done.add(outputbase)
    return already_done


def check_domtable_files(run, cluster_base_names, cached_domtable_files):
    # domtableFiles: all domtable files corresponding to the input files
    existing_domtable_files = set()
    for name in cluster_base_names:
        existing_domtable_files.add(os.path.join(run.directories.domtable, name+".domtable"))

    # domtableBases: the actual set of input files with coresponding domtable files
    domtable_bases = cached_domtable_files.intersection(existing_domtable_files)

    # Verify that all input files have a corresponding domtable file
    if len(existing_domtable_files - domtable_bases) > 0:
        skipped_files = existing_domtable_files - domtable_bases
        logging.error("The following files did NOT have their domains predicted: ")
        for skipped_file in skipped_files:
            logging.error(skipped_file)
        sys.exit(1)


def log_unprocessed_domtable_files(files_to_process, first_run):
    if first_run:
        logging.info(" Processing %d domtable files", len(files_to_process))
        return
    else:
        if len(files_to_process) == 0:
            logging.info(" All domtable files had already been processed")
        elif len(files_to_process) < 20:
            logging.info(" The following domtable files had not been processed:")
            for unprocessed_domtable_file in files_to_process:
                logging.info(unprocessed_domtable_file.split(os.sep)[-1].split('.')[:-1][0])
        else:
            logging.info(" %d domtable files will be processed", len(files_to_process))

def get_cached_pfd_files(run):
    file_path = os.path.join(run.directories.pfd, "*.pfd")
    pfd_files = glob(file_path)
    return set(pfd_files)


def check_pfd_files(run, cluster_base_names):
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
