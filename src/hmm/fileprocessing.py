import os
import sys
from glob import glob

def find_unprocessed_files(run, fasta_files):
    ALREADY_DONE = set()
    for fasta in fasta_files:
        outputbase = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
        outputfile = os.path.join(run.directories.domtable, outputbase + '.domtable')
        if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
            # verify domtable content
            with open(outputfile, "r") as domtablefile:
                for line in domtablefile.readlines():
                    if line.startswith("# Option settings:"):
                        linecols = line.split()
                        if "hmmscan" in linecols and "--domtblout" in linecols:
                            ALREADY_DONE.add(fasta)
                            break

    TASK_SET = fasta_files - ALREADY_DONE
    if len(TASK_SET) == 0:
        print(" All fasta files had already been processed")
    elif len(ALREADY_DONE) > 0:
        if len(TASK_SET) < 20:
            TASKS = [x.split(os.sep)[-1].split(".")[:-1] for x in TASK_SET]
            print(" Warning! The following NEW fasta file(s) will be processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in TASK_SET)))
        else:
            print(" Warning: {} NEW fasta files will be processed".format(len(TASK_SET)))
    else:
        print(" Predicting domains for {} fasta files".format(str(len(fasta_files))))
    
    return TASK_SET


def verify_hmm_fasta(run, base_names):
    # All available fasta files (could be more than it should if reusing output folder)
    all_fasta_files = set(glob(os.path.join(run.directories.bgc_fasta, "*.fasta")))

    # fastaFiles: all the fasta files that should be there
    # (i.e. correspond to the input files)
    fasta_files = set()
    for name in base_names:
        fasta_files.add(os.path.join(run.directories.bgc_fasta, name+".fasta"))

    # fastaBases: the actual fasta files we have that correspond to the input
    fasta_bases = all_fasta_files.intersection(fasta_files)

    # Verify that all input files had their fasta sequences extracted
    if len(fasta_files - fasta_bases) > 0:
        print("Error! The following files did NOT have their fasta sequences extracted: ")
        unextracted_files = fasta_files - fasta_bases
        for unextracted_file in unextracted_files:
            print(unextracted_file)
        sys.exit()
    
    return all_fasta_files