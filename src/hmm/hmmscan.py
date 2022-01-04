import os
import sys
import subprocess
import multiprocessing

from src.hmm.io import domtable_parser
from src.pfam.io import write_pfd
from src.pfam.misc import check_overlap

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

def run_hmmscan_multi_threaded(RUN, TASK_SET):
    POOL = multiprocessing.Pool(RUN.options.cores, maxtasksperchild=1)
    for fastaFile in TASK_SET:
        task_args = (fastaFile, RUN.directories.pfam, RUN.directories.domtable, RUN.options.verbose)
        POOL.apply_async(runHmmScan, args=task_args)
    POOL.close()
    POOL.join()

def runHmmScan(fastaPath, hmmPath, outputdir, verbose):
    """ Runs hmmscan command on a fasta file with a single core to generate a
    domtable file"""
    hmmFile = os.path.join(hmmPath,"Pfam-A.hmm")
    if os.path.isfile(fastaPath):
        name = ".".join(fastaPath.split(os.sep)[-1].split(".")[:-1])
        outputName = os.path.join(outputdir, name+".domtable")
        
        hmmscan_cmd = "hmmscan --cpu 0 --domtblout {} --cut_tc {} {}".format(outputName, hmmFile, fastaPath)
        if verbose == True:
            print("   " + hmmscan_cmd)
        subprocess.check_output(hmmscan_cmd, shell=True)

    else:
        sys.exit("Error running hmmscan: Fasta file " + fastaPath + " doesn't exist")


def parseHmmScan(hmmscanResults, pfd_folder, pfs_folder, overlapCutoff, verbose, genbankDict, clusters, baseNames, gbk_files, sampleDict, mibig_set):
    outputbase = ".".join(hmmscanResults.split(os.sep)[-1].split(".")[:-1])
    # try to read the domtable file to find out if this gbk has domains. Domains
    # need to be parsed into fastas anyway.
    if os.path.isfile(hmmscanResults):
        pfd_matrix = domtable_parser(outputbase, hmmscanResults)
        
        # get number of domains to decide if this BGC should be removed
        num_domains = len(pfd_matrix)

        if num_domains > 0:
            if verbose:
                print("  Processing domtable file: " + outputbase)

            # check_overlap also sorts the filtered_matrix results and removes
            # overlapping domains, keeping the highest scoring one
            filtered_matrix, domains = check_overlap(pfd_matrix,overlapCutoff)
            
            # Save list of domains per BGC
            pfsoutput = os.path.join(pfs_folder, outputbase + ".pfs")
            with open(pfsoutput, 'w') as pfs_handle:
                pfs_handle.write(" ".join(domains))
            
            # Save more complete information of each domain per BGC
            pfdoutput = os.path.join(pfd_folder, outputbase + ".pfd")
            with open(pfdoutput,'w') as pfd_handle:
                write_pfd(pfd_handle, filtered_matrix)
        else:
            # there aren't any domains in this BGC
            # delete from all data structures
            print("  No domains where found in {}.domtable. Removing it from further analysis".format(outputbase))
            info = genbankDict.get(outputbase)
            clusters.remove(outputbase)
            baseNames.remove(outputbase)
            gbk_files.remove(info[0])
            for sample in info[1]:
                sampleDict[sample].remove(outputbase)
            del genbankDict[outputbase]
            if outputbase in mibig_set:
                mibig_set.remove(outputbase)
            
    else:
        sys.exit("Error: hmmscan file " + outputbase + " was not found! (parseHmmScan)")

    return("")
