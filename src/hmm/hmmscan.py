import os
import sys
import subprocess


from src.hmm.io import domtable_parser
from src.pfam.io import write_pfd
from src.pfam.misc import check_overlap

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
