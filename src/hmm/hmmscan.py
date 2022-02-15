import logging
import os
import sys
import subprocess
import multiprocessing

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
        logging.error("  Could not find file %s", dom_file)
        sys.exit(1)
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
                    logging.warning("No gene ID in %s", gbk)
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

def write_pfd(pfd_handle, matrix):
    for row in matrix:
        row = "\t".join(row)
        pfd_handle.write(row+"\n")

    pfd_handle.close()


def run_hmmscan_async(run, task_set):
    pool = multiprocessing.Pool(run.options.cores, maxtasksperchild=1)
    for fasta_file in task_set:
        task_args = (fasta_file, run.directories.pfam, run.directories.domtable, run.options.verbose)
        pool.apply_async(run_hmmscan, args=task_args)
    pool.close()
    pool.join()

def run_hmmscan(fasta_path, hmm_path, outputdir, verbose):
    """ Runs hmmscan command on a fasta file with a single core to generate a
    domtable file"""
    hmm_file = os.path.join(hmm_path, "Pfam-A.hmm")
    if os.path.isfile(fasta_path):
        name = ".".join(fasta_path.split(os.sep)[-1].split(".")[:-1])
        output_name = os.path.join(outputdir, name+".domtable")

        hmmscan_cmd = "hmmscan --cpu 0 --domtblout {} --cut_tc {} {}".format(output_name, hmm_file, fasta_path)

        logging.debug("   %s", hmmscan_cmd)
        subprocess.check_output(hmmscan_cmd, shell=True)

    else:
        logging.error("Error running hmmscan: Fasta file %s doesn't exist", fasta_path)
        sys.exit(1)

def parse_hmmscan(hmm_scan_results, pfd_folder, pfs_folder, overlap_cutoff, verbose, genbank_dict, clusters, base_names, mibig_set):
    sample_dict = {} # {sampleName:set(bgc1,bgc2,...)}
    gbk_files = [] # raw list of gbk file locations
    for (cluster, (path, cluster_sample)) in genbank_dict.items():
        gbk_files.append(path)
        for sample in cluster_sample:
            clusters_in_sample = sample_dict.get(sample, set())
            clusters_in_sample.add(cluster)
            sample_dict[sample] = clusters_in_sample

    outputbase = ".".join(hmm_scan_results.split(os.sep)[-1].split(".")[:-1])
    # try to read the domtable file to find out if this gbk has domains. Domains
    # need to be parsed into fastas anyway.
    if os.path.isfile(hmm_scan_results):
        pfd_matrix = parse_domtable(outputbase, hmm_scan_results)

        # get number of domains to decide if this BGC should be removed
        num_domains = len(pfd_matrix)

        if num_domains > 0:
            logging.debug("  Processing domtable file: %s", outputbase)

            # check_overlap also sorts the filtered_matrix results and removes
            # overlapping domains, keeping the highest scoring one
            filtered_matrix, domains = check_overlap(pfd_matrix, overlap_cutoff)

            # Save list of domains per BGC
            pfsoutput = os.path.join(pfs_folder, outputbase + ".pfs")
            with open(pfsoutput, 'w') as pfs_handle:
                pfs_handle.write(" ".join(domains))

            # Save more complete information of each domain per BGC
            pfdoutput = os.path.join(pfd_folder, outputbase + ".pfd")
            with open(pfdoutput, 'w') as pfd_handle:
                write_pfd(pfd_handle, filtered_matrix)
        else:
            # there aren't any domains in this BGC
            # delete from all data structures
            logging.info("  No domains were found in %s.domtable. Removing it from further analysis", outputbase)
            info = genbank_dict.get(outputbase)
            clusters.remove(outputbase)
            base_names.remove(outputbase)
            gbk_files.remove(info[0])
            for sample in info[1]:
                sample_dict[sample].remove(outputbase)
            del genbank_dict[outputbase]
            if outputbase in mibig_set:
                mibig_set.remove(outputbase)

    else:
        logging.error("hmmscan file %s was not found! (parseHmmScan)", outputbase)
        sys.exit(1)
