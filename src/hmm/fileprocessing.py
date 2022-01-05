import os
import sys
from glob import glob

def get_cached_fasta_files(run):
    # All available fasta files (could be more than it should if reusing output folder)
    file_path = os.path.join(run.directories.bgc_fasta, "*.fasta")
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
        print("Error! The following files did NOT have their fasta sequences extracted: ")
        unextracted_files = existing_fasta_files - fasta_bases
        for unextracted_file in unextracted_files:
            print(unextracted_file)
        sys.exit()

def get_fasta_files_to_process(run, fasta_files):
    if run.options.force_hmmscan:
        print(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
        return fasta_files

    already_done = set()
    for fasta in fasta_files:
        outputbase = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
        outputfile = os.path.join(run.directories.domtable, outputbase + '.domtable')

        # ignore directories or empty files
        if not os.path.isfile(outputfile) or os.path.getsize(outputfile) == 0:
            continue

        # verify domtable content
        with open(outputfile, "r") as domtablefile:
            for line in domtablefile.readlines():
                if line.startswith("# Option settings:"):
                    linecols = line.split()
                    if "hmmscan" in linecols and "--domtblout" in linecols:
                        already_done.add(fasta)
                        break

    task_set = fasta_files - already_done
    if len(task_set) == 0:
        print(" All fasta files had already been processed")
    elif len(already_done) > 0:
        if len(task_set) < 20:
            tasks = ", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in task_set)
            print(" Warning! The following NEW fasta file(s) will be processed: {}".format(tasks))
        else:
            print(" Warning: {} NEW fasta files will be processed".format(len(task_set)))
    else:
        print(" Predicting domains for {} fasta files".format(str(len(fasta_files))))

    return task_set

def get_cached_domtable_files(run):
    file_path = os.path.join(run.directories.domtable, "*.domtable")
    domtable_files = glob(file_path)
    return set(domtable_files)

def get_domtable_files_to_process(run, domtable_files):
    # return early if force_hmmscan is true
    if run.options.force_hmmscan:
        return domtable_files

    already_done = get_processed_domtable_files(run, domtable_files)

    # the domtable files to be processed is the difference between all, and those which were 
    # already done
    unprocessed_domtable_files = domtable_files - already_done
    first_run = len(domtable_files) == 0

    log_unprocessed_domtable_files(unprocessed_domtable_files, first_run)

    return unprocessed_domtable_files

def get_processed_domtable_files(run, domtable_files):
    # otherwise find what still needs to be done
    already_done = set()
    for domtable in domtable_files:
        outputbase = ".".join(domtable.split(os.sep)[-1].split(".")[:-1])
        outputfile = os.path.join(run.directories.pfd, outputbase + '.pfd')
        if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
            already_done.add(domtable)
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
        print("Error! The following files did NOT have their domains predicted: ")
        for skipped_file in skipped_files:
            print(skipped_file)
        sys.exit()


def log_unprocessed_domtable_files(files_to_process, first_run):
    if first_run:
        print(" Processing {} domtable files".format(str(len(files_to_process))))
        return
    else:
        if len(files_to_process) == 0:
            print(" All domtable files had already been processed")
        elif len(files_to_process) < 20:
            print(" Warning! The following domtable files had not been processed:")
            for unprocessed_domtable_file in files_to_process:
                print(unprocessed_domtable_file.split(os.sep)[-1].split('.')[:-1][0])
        else:
            print(" Warning: {} domtable files will be processed".format(str(len(files_to_process))))
    
