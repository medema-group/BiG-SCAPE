import os
import sys
import subprocess
import multiprocessing

from src.pfam.io import write_pfd
from src.pfam.misc import check_overlap

def parse_domtable(gbk, dom_file):
    """Parses the domain table output files from hmmscan"""

##example from domain table output:
    # target name        accession   tlen query name                                    accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
    #------------------- ---------- -----                          -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
    #Lycopene_cycl        PF05834.8    378 loc:[0:960](-):gid::pid::loc_tag:['ctg363_1'] -            320   3.1e-38  131.7   0.0   1   1   1.1e-40   1.8e-36  126.0   0.0     7   285    33   295    31   312 0.87 Lycopene cyclase protein

    #pfd_matrix columns: gbk filename - score - gene id - first coordinate - second coordinate - pfam id - domain name -start coordinate of gene - end coordinate of gene - cds header

    pfd_matrix = []
    try:
        dom_handle = open(dom_file, 'r')
    except IOError:
        print("  Error! Could not find file " + dom_file)
    else:
        for line in dom_handle:
            if line[0] != "#":
                splitline = line.split()
                pfd_row = []
                pfd_row.append(gbk)         #add clustername or gbk filename

                pfd_row.append(splitline[13]) #add the score

                header_list = splitline[3].split(":")
                try:
                    pfd_row.append(header_list[header_list.index("gid")+1]) #add gene ID if known
                except ValueError:
                    print("No gene ID in " + gbk)
                    pfd_row.append('')

                pfd_row.append(splitline[19])#first coordinate, env coord from
                pfd_row.append(splitline[20])#second coordinate, env coord to

                #===================================================================
                # loc_split = header_list[2].split("]") #second coordinate (of CDS) and the direction
                # pfd_row.append(loc_split[1]) #add direction
                #===================================================================

                pfd_row.append(splitline[1]) #pfam id
                pfd_row.append(splitline[0]) #domain name
                pfd_row.append(header_list[header_list.index("loc")+1]) #start coordinate of gene
                pfd_row.append(header_list[header_list.index("loc")+2]) #end coordinate of gene

                pfd_row.append(splitline[3])#cds header
                pfd_matrix.append(pfd_row)

    return pfd_matrix


def run_hmmscan_async(RUN, TASK_SET):
    POOL = multiprocessing.Pool(RUN.options.cores, maxtasksperchild=1)
    for fastaFile in TASK_SET:
        task_args = (fastaFile, RUN.directories.pfam, RUN.directories.domtable, RUN.options.verbose)
        POOL.apply_async(run_hmmscan, args=task_args)
    POOL.close()
    POOL.join()

def run_hmmscan(fastaPath, hmmPath, outputdir, verbose):
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

def parse_hmmscan_async():
    # copied from bigscape.py, unused
    # TODO: evaluate if still needed, reimplement?
    # If using the multiprocessing version and outputbase doesn't have any
    #  predicted domains, it's not as easy to remove if from the analysis
    #  (probably because parseHmmScan only has a copy of clusters et al?)
    # Using serialized version for now. Probably doesn't have too bad an impact
    #pool = Pool(cores,maxtasksperchild=32)
        #task_args = (domtableFile,output_folder,options.domain_overlap_cutoff)
        #pool.apply_async(parseHmmScan, args = task_args)
    #pool.close()
    #pool.join()
    return

def parse_hmmscan(hmmscanResults, pfd_folder, pfs_folder, overlapCutoff, verbose, genbankDict, clusters, baseNames, mibig_set):
    sampleDict = {} # {sampleName:set(bgc1,bgc2,...)}
    gbk_files = [] # raw list of gbk file locations
    for (cluster, (path, clusterSample)) in genbankDict.items():
        gbk_files.append(path)
        for sample in clusterSample:
            clustersInSample = sampleDict.get(sample, set())
            clustersInSample.add(cluster)
            sampleDict[sample] = clustersInSample

    outputbase = ".".join(hmmscanResults.split(os.sep)[-1].split(".")[:-1])
    # try to read the domtable file to find out if this gbk has domains. Domains
    # need to be parsed into fastas anyway.
    if os.path.isfile(hmmscanResults):
        pfd_matrix = parse_domtable(outputbase, hmmscanResults)

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
