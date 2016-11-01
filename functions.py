#!/usr/bin/env python

"""
BiG-SCAPE

PI: Marnix Medema

Developers:
Marley Yeong                    marleyyeong@live.nl
Jorge Navarro
Emmanuel (Emzo) de los Santos

Functions used by bigscape.py

"""

import os
import subprocess
import sys
from encodings import gbk
from Bio import SeqIO
import copy
import math

global verbose
verbose = False

def frange(start, stop, step):
    """Range function with float increments"""
    
    i = start
    while i <= stop:
        yield i
        i += step


def get_anchor_domains(filename):
    """Get the anchor/marker domains from a txt file.
    This text file should contain one Pfam id per line"""

    domains = []
    
    try:
        handle = open(filename, 'r')
    except IOError:
        print "You have not provided the anchor_domains.txt file."
        print "if you want to make use of the anchor domains in the DDS distance metric,\
        make a file that contains a Pfam domain on each line."
        return []
        
    raw_domains = handle.readlines()
    for line in raw_domains:
        domains.append(line.strip())
    handle.close()
    return domains
        

def get_domain_list(filename):
    """Convert the Pfam string in the .pfs files to a Pfam list"""
    handle = open(filename, 'r')
    domains_string = handle.readline().strip()
    domains = domains_string.split(" ")
    handle.close()
    
    return domains


def get_all_features_of_type(gb_file, types):
    "Return all features of the specified types for a seq_record"
    
    handle = open(gb_file, "r")
    seq_record = SeqIO.read(handle, "genbank")
    
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    handle.close()
            
    return features


def get_hmm_output_files(output_folder):
    """Finds all .domtable files in the output folder and its child folders"""
    
    hmm_table_list = []
    for dirpath, dirnames, filenames in os.walk(str(output_folder) + "/"):
        for fname in filenames:
            if fname.split(".")[-1] == "domtable":
                try:
                    if open(output_folder + "/" + fname, "r").readlines()[3][0] != "#": #if this is false, hmmscan has not found any domains in the sequence
                        hmm_table_list.append(fname)
                        if verbose == True:
                            print fname
                    else:
                        print output_folder + "/" + fname, "seems to be an empty domtable file"
                    
                except IndexError:
                    print "IndexError on file", output_folder + "/" + fname

    return hmm_table_list


def check_overlap(pfd_matrix, overlap_cutoff):
    """Check if domains overlap for a certain overlap_cutoff.
     If so, remove the domain(s) with the lower score."""
    
    delete_list = []
    for i in range(len(pfd_matrix)-1):
        for j in range(i+1, len(pfd_matrix)):
            row1 = pfd_matrix[i]
            row2 = pfd_matrix[j]
            
            #check if we are the same CDS
            if row1[-1] == row2[-1]:
                #check if there is overlap between the domains
                if no_overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4])) == False:
                    overlapping_nucleotides = overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4]))
                    overlap_perc_loc1 = overlap_perc(overlapping_nucleotides, int(row1[4])-int(row1[3]))
                    overlap_perc_loc2 = overlap_perc(overlapping_nucleotides, int(row2[4])-int(row2[3]))
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[1]) > float(row2[1]): #see which has a better score
                            delete_list.append(row2)
                        elif float(row1[1]) < float(row2[1]):
                            delete_list.append(row1)

    for lst in delete_list:
        try:
            pfd_matrix.remove(lst)
        except ValueError:
            pass
        
    # for some reason, some coordinates in genbank files have ambiguous 
    # starting/ending positions for the CDS. In this case we need the 
    # absolute start position of each domain so the loci start coordinate
    # needs to be checked
    absolute_start_positions = [] # we have to take care of strandedness of the origin gene
    for row in pfd_matrix:
        if row[7][0] == "<" or row[7][0] == ">":
            row[7] = row[7][1:]
        if row[8][0] == "<" or row[8][0] == ">":
            row[8] = row[8][1:]
        
        loci_start = int(row[7])
        loci_end = int(row[8])
        domain_start = int(row[3])
        domain_end = int(row[4])
        width = 3*(domain_end - domain_start)
        strand = row[9].split(":")[-1]
        if strand == "+":
            absolute_start_positions.append(3*domain_start + loci_start)
        else:
            absolute_start_positions.append(loci_end - 3*domain_start - width)
        
    try:
        # sorted by absolute position: env coordinate start + gene locus start
        pfd_matrix = [x for (y,x) in sorted(zip(absolute_start_positions,pfd_matrix), key=lambda pair:pair[0])]
    except ValueError as e:
        print("Something went wrong during the sorting of the fourth column (ValueError): " + str(e))
        print(pfd_matrix)
        sys.exit()
        
    domains = []
    for row in pfd_matrix:
        domains.append(row[5]) #save the pfam domains for the .pfs file

    return pfd_matrix, domains                       
                        

def write_pfs(pfs_handle, domains):
    for domain in domains:
        pfs_handle.write(domain+" ")
    pfs_handle.close()
    

def write_pfd(pfd_handle, matrix):
    for row in matrix:
        row = "\t".join(row)
        pfd_handle.write(row+"\n")
        
    pfd_handle.close() 
    

def dct_writeout(handle, dct):
    for key in dct.keys():
        handle.write(key+"\n"+dct[key]+"\n")
    handle.close()


def check_data_integrity(gbk_files):
    """Perform some integrity checks on the input gbk files."""

    discarded_set = set()

    # check for potential errors while reading the gbk files
    # Structure must remain intact after removing offending files
    # -> The following used to be necessary before storing sample information
    # in genbankDict. It could eventually be removed/simplified
    gbk_files_check = copy.deepcopy(gbk_files)
    samples_for_deletion = []
    for sample in range(len(gbk_files_check)):
        for f in gbk_files_check[sample]:
            handle = open(f, "r")
            try:
                SeqIO.read(handle, "genbank")
            except ValueError as e:
                print("   Error with file " + f + ": \n    '" + str(e) + "'")
                print("    (This file will be excluded from the analysis)")
                discarded_set.add(f)
                if len(gbk_files[sample]) > 1:
                    gbk_files[sample].remove(f)
                else:
                    # this is a one-element list, mark it for late removal
                    samples_for_deletion.append(sample)
                pass
            handle.close()
    if len(samples_for_deletion) > 0:
        for s in sorted(samples_for_deletion, reverse=True):
            del gbk_files[s]
            
    # check for duplication
    duplication = False
    del gbk_files_check[:]
    gbk_files_check = list(iterFlatten(gbk_files))
    for file in gbk_files_check:
        name = file.split(os.sep)[-1]
        file_occ = 0
        for cfile in gbk_files_check:
            
            cname = cfile.split(os.sep)[-1]
            if name == cname:
                file_occ += 1
                if file_occ > 1:
                    duplication = True
                    print "duplicated file at:", cfile 
    
    # Without proper checking downstream, allowing duplicate files results in
    # problems: domain's sequences could be duplicated and the length's mismatch
    # raises an exception in sequence similarity scoring. 
    # This probably is not an issue anymore with the new genbankDict structure
    if duplication == True:
        print "There was duplication in the input files, if this is not intended remove them."
        cont = raw_input("Continue anyway? Y/N ")
        if cont.lower() == "n":
            sys.exit()
            
    # The possibility of not having files to analyze was checked for in get_gbk_files
    # so if the list is empty now, it's due to issues in the gbk files
    if len(gbk_files_check) < 2:
        sys.exit("Due to errors in the input files, there are no files left for the analysis (" + str(len(gbk_files_check)) + ")")
       
    return discarded_set
            

def hmmscan(pfam_dir, fastafile, outputdir, name, cores):
    """Runs hmmscan"""
    #removed --noali par

    hmmscan_cmd = "hmmscan --cpu " + str(cores) + " --domtblout " + os.path.join(outputdir, name+".domtable") + " --cut_tc " + os.path.join(pfam_dir,"Pfam-A.hmm") + " " + str(fastafile)
    if verbose == True:
        print("\t"+hmmscan_cmd)
    
    subprocess.check_output(hmmscan_cmd, shell=True)
    
    
def get_domains(filename):
    handle = open(filename, 'r')
    domains = []
    for line in handle:
        if line[0] != "#":
            domains.append(filter(None, line.split(" "))[1])
            
    return domains


def no_overlap(locA1, locA2, locB1, locB2):
    """Return True if there is no overlap between two regions"""
    if locA1 < locB1 and locA2 < locB1:
        return True
    elif locA1 > locB2 and locA2 > locB2:
        return True
    else:
        return False


def overlap_perc(overlap, len_seq):
    return float(overlap) / len_seq
    

def overlap(locA1, locA2, locB1, locB2):
    """Returns the amount of overlapping nucleotides"""

    if locA1 < locB1:
        cor1 = locA1
    else:
        cor1 = locB1

    if locA2 > locB2:
        cor2 = locA2
    else:
        cor2 = locB2

    total_region = cor2 - cor1
    sum_len = (locA2 - locA1) + (locB2 - locB1)

    return sum_len - total_region

  
def calc_perc_identity(seq1, seq2, spec_domain, spec_domain_nest, domain):
    """Percent Identity = Matches/Length of aligned region (with gaps) (sequence similarity!)
    Positions where both characters are gaps (dashes) are not counted in the Length"""

    length = 0
    matches = 0
    
    #Sequences *should* have the same length because they come from an MSA
    # but if sequences were incorrectly appended more than once to the same 
    # domain sequence, there would be problems here
    if len(seq1) != len(seq2):
        print("\tWARNING: mismatch in sequences' lengths while calculating sequence identity") 
        print("\t Domain: " + domain)
        print("\t  Specific domain 1: " + spec_domain + " len: " + str(len(seq1)))
        print("\t  Specific domain 2: " + spec_domain_nest + " len: " + str(len(seq2)))
        print("\t trying to continue with shortest length...")
        
    seq_length = min(len(seq1), len(seq2))
    
    for pos in range(seq_length): 
        if seq1[pos] == seq2[pos]:
            if seq1[pos] != "-":
                matches += 1
                length += 1
        else:
            length += 1

    return float(matches) / float(length), length    


def BGC_dic_gen(filtered_matrix):
    """Generates the: { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } part of the BGCs variable."""
    bgc_dict = {}
    for row in filtered_matrix:
        try: #Should be faster than performing if key in dictionary.keys()
            bgc_dict[row[5]]
            bgc_dict[row[5]].append(str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4]))
        except KeyError: #In case of this error, this is the first occurrence of this domain in the cluster
           bgc_dict[row[5]]=[str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4])]
            
    return bgc_dict

    
def save_domain_seqs(filtered_matrix, fasta_dict, domains_folder, output_folder, outputbase):
    """Write fasta sequences for the domains in the right pfam-domain file"""
    for row in filtered_matrix:
        domain = row[5]
        seq = fasta_dict[">"+str(row[-1].strip())] #access the sequence by using the header
        

        domain_file = open(os.path.join(output_folder, domains_folder, domain + ".fasta"), 'a') #append to existing file
        domain_file.write(">" + str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4]) \
        + "\n" + str(seq)[int(row[3]):int(row[4])] + "\n") #only use the range of the pfam domain within the sequence
        
#===============================================================================
# 
#         if str(seq)[int(row[3]):int(row[4])] == "":
#             print seq
#             print row[3], row[4]
#             print str(row[-1].strip())
#             print outputbase
#===============================================================================
            
        domain_file.close()


def network_parser(network_file, Jaccardw, DDSw, GKw):
    network = {}
    
    try:
        with open(network_file, "r") as handle:
            next(handle) # ignore first line with column names
            for line in handle:
                if line[0] != "#":
                    strippedline = line.strip().split("\t")
                    
                    # don't read the combined and shared group (last two columns)
                    # as these are generated by write_network_matrix()
                    network[strippedline[0], strippedline[1]] = strippedline[2:6]
                    network[strippedline[0], strippedline[1]].extend(float(x) for x in strippedline[6:-2]) 
        handle.close()
    except IOError:
        sys.exit("Error: Cannot open file " + network_file)
    
    # re-calculate raw distance with potentially new weights
    for (a, b) in network:
        Jaccard = network[a,b][7]
        DDS = network[a,b][8]
        GK = network[a,b][9]
        distance = 1- (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK)
        if distance <= 0:
            logscore = float("inf")
        else:
            logscore = -log(distance, 2)
        sqrd_similarity = (1-distance)**2
        
        network[a,b][4] = logscore
        network[a,b][5] = distance
        network[a,b][6] = sqrd_similarity
        
    return network


def fasta_parser(handle):
    """Parses a fasta file, and stores it in a dictionary.
    Only works if there are no duplicate fasta headers in the fasta file"""
    
    fasta_dict = {}
    header = ""
    for line in handle:
        if line[0] == ">":
            header=line.strip()
        else:
            try:
                fasta_dict[header] += line.strip()
            except KeyError:
                fasta_dict[header] = line.strip()

    return fasta_dict


def get_domain_fastas(domain_folder, output_folder):
    """Finds the pfam domain fasta files"""
    domain_fastas = []
    dirpath_ = ""
    for dirpath, dirnames, filenames in os.walk(os.path.join(output_folder, domain_folder)):
        for fname in filenames:
            if ".fasta" in fname and "hat2" not in fname:
                domain_fastas.append(os.path.join(dirpath, fname))
                if verbose == True:
                    print fname
                    
    return domain_fastas


def iterFlatten(root):
    """'Flattens' a matrix by collapsing the lists into one list"""
    if isinstance(root, (list, tuple)):
        for element in root:
            for e in iterFlatten(element):
                yield e
    else:
        yield root


def distout_parser(distout_file):
    """returns similarity values, for domains in the following format  { ('specific_domain_name_1',
    'specific_domain_name_2'): (sequence_identity, alignment_length), ... }"""
    
    try:
        hat2_handle = open(distout_file, 'r')
    except IOError:
        return {}
    
    domain_pairs_dict = {}    
    linecounter = 0
    seqsdict = {}
    distances = [] #will be of length numberof_seqs * (numberof_seqs-1) / 2
    numberof_seqs = 0
    for line in hat2_handle:
        linecounter += 1
        if linecounter == 2: #contains number of sequences
            numberof_seqs = int(line.replace(" ", "").strip())
            
        elif linecounter >= 4 and linecounter <= 3 + numberof_seqs:
            try:
                #seq_number = int(re.search(r' \d*\. ', str(line.split("=")[0])).group(0).replace(".", "").replace(" ", ""))
                seq_number = int(line.split("=")[0].replace(" ", "").replace(".", ""))
            except AttributeError:
                print "something went wrong during the import of distout file: ", str(distout_file)
                
            
            seqsdict[seq_number] = "".join(line.split("=")[1:]).strip()#in case the header contains an = sign

        elif linecounter > 3 + numberof_seqs:
            distances += line.strip().split(" ")

    keys=[]
    if len(distances) != (numberof_seqs * (numberof_seqs-1)) / 2.0:
        print "something went horribly wrong in importing the distance matrix"
    else:
        print "distance matrix imported correctly"
        keys = seqsdict.keys()

    keys_queue = []
    for key in keys:
        keys_queue.append(key)

    tuples = []
    for key in keys:
        keys_queue.remove(key)
        for key_queue in keys_queue:
            tuples.append((key, key_queue))
            
    for tupl in range(len(tuples)):
        ##    { ('specific_domain_name_1',
        ##    'specific_domain_name_2'): (sequence_identity, alignment_length), ... }
        #1-distance is a representation of the sequence_identity
        domain_pairs_dict[tuple(sorted([seqsdict[tuples[tupl][0]], seqsdict[tuples[tupl][1]]]))] = (1-float(distances[tupl]), 0)

    return domain_pairs_dict


def write_network_matrix(matrix, cutoff, filename, include_disc_nodes):
    networkfile = open(filename, 'w')
    clusters = [] # will contain the names of clusters that have an edge value lower than the threshold
    networkfile.write("clustername1\tclustername2\tgroup1\tdefinition\tgroup2\tdefinition\t-log2score\traw distance\tsquared similarity\tJaccard index\tDDS index\tGK index\tcombined group\tshared group\n")
    
    for (gc1, gc2) in matrix.keys():
        row = [gc1, gc2]
        for i in matrix[gc1, gc2]:
            row.append(i)

        temprow = []
        #if both clusters have the same group, this group will be annotated in the last column of the network file
        for i in row:
            temprow.append(i)
        
        if float(temprow[7]) <= float(cutoff):
            clusters.append(row[0])
            clusters.append(row[1])
            
            if row[2] != "" and row[4] != "": #group1, group2
                temprow.append(" - ".join(sorted([str(row[2]),str(row[4])])))
            elif row[4] != "":
                temprow.append(str(row[4]))
            elif row[2] != "":
                temprow.append(str(row[2]))
            else:
                 temprow.append(str("NA"))
            
            if row[2] == row[4]:
                temprow.append(row[2])
            else:
                temprow.append("")
                
            networkfile.write("\t".join(map(str,temprow)) + "\n")

    # matrix[gc1, gc2] =
    # row:   0      1    2    3       4       5      6      7    8   9   [   10      11   <- these two are written directly
    #       grp1  def1 grp2  def2  -logScr  rawD  sqrtSim  Jac  DDS  GK  [combGrp  ShrdGrp
    if include_disc_nodes == True:  
        #Add the nodes without any edges, give them an edge to themselves with a distance of 0 
        clusters = set(clusters)
        passed_clusters = []
        for (gc1, gc2) in matrix.keys():
            if gc1 not in clusters and gc1 not in passed_clusters:
                networkfile.write("\t".join([gc1, gc1, matrix[gc1, gc2][0], matrix[gc1, gc2][1], "0", "0", "1", "0", "0", "0", "", ""]) + "\n")
                passed_clusters.append(gc1)
            
            if gc2 not in clusters and gc2 not in passed_clusters:
                networkfile.write("\t".join([gc2, gc2, matrix[gc1, gc2][2], matrix[gc1, gc2][3], "0", "0", "1", "0", "0", "0", "", ""]) + "\n")
                passed_clusters.append(gc2)
            
    networkfile.close()
                    

def domtable_parser(gbk, dom_file):
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
                splitline = filter(None, line.split(" "))
                pfd_row = []
                pfd_row.append(gbk)         #add clustername or gbk filename
                
                pfd_row.append(splitline[13]) #add the score

                header_list = splitline[3].split(":")
                try:
                    pfd_row.append(header_list[header_list.index("gid")+1]) #add gene ID if known
                except ValueError:
                    print "No gene ID in ", gbk
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


def hmm_table_parser(gbk, hmm_table):
##AB050629.gbk    score   yxjC    1    1167    +    PF03600    CitMHS

##example from hmm table output:
##Thiolase_N           PF00108.19 loc:[2341:3538](+):gid::pid::loc_tag:['ctg4508_7'] -
##8.2e-90  300.5   2.3   1.2e-89  300.0   2.3   1.2   1   0   0   1   1   1   1 Thiolase, N-terminal domain

    pfd_matrix = []

    handle = open(hmm_table, 'r')
    for line in handle:
        
        if line[0] != "#":
            
            splitline = filter(None, line.split(" "))
            pfd_row = []
            pfd_row.append(gbk)
            pfd_row.append(splitline[5]) #add the score

##example of header_list ['loc', '[2341', '3538](+)', 'gid', '', 'pid', '', 'loc_tag', "['ctg4508_7"]]

            header_list = splitline[2].split(":")
            pfd_row.append(header_list[header_list.index("gid")+1])
            pfd_row.append(header_list[1].replace("[", "")) #first coordinate
            loc_split = header_list[2].split("]") #second coordinate and the direction
            pfd_row.append(loc_split[0])
            pfd_row.append(loc_split[1])

            pfd_row.append(splitline[1])
            pfd_row.append(splitline[0])
            
            pfd_matrix.append(pfd_row)

    return pfd_matrix


def write_parameters(output_folder, options):
    """Write a file with all the details of the run.
    Will overwrite previous versions"""
    
    pf = open(os.path.join(output_folder,"parameters.txt"), "w")
    
    pf.write("Input directory:\t" + options.inputdir + "\n")
    pf.write("Cores:\t" + str(options.cores) + "\n")
    
    pf.write("Verbose:\t" + ("True" if options.verbose else "False"))
    if not options.verbose:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("Include disc nodes:\t" + ("True" if options.include_disc_nodes else "False"))
    if not options.include_disc_nodes:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("Domain overlap cutoff:\t" + str(options.domain_overlap_cutoff))
    if options.domain_overlap_cutoff == 0.1:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Minimum size of BGC:\t" + str(options.min_bgc_size))
    if options.min_bgc_size == 0:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    pf.write("Sequence distance networks:\t\"" + str(options.seqdist_networks) + "\"")
    if options.seqdist_networks == "A":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("Domain distance networks:\t\"" + str(options.domaindist_networks) + "\"")
    if options.domaindist_networks == "":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    pf.write("Jaccard weight:\t" + str(options.Jaccardw))
    if options.Jaccardw == 0.2:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("DDS weight:\t" + str(options.DDSw))
    if options.DDSw == 0.75:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("GK weight:\t" + str(options.GKw))
    if options.GKw == 0.05:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Anchor domain weight:\t" + str(options.anchorweight))
    if options.anchorweight == 0.1:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    #pf.write("Output folder for domain fasta files:\t" + options.domainsout)
    #if options.domainsout == "domains":
        #pf.write("\t(default)\n")
    #else:
        #pf.write("\n")
        
    pf.write("Location of Pfam files:\t" + options.pfam_dir)
    if options.pfam_dir == os.path.dirname(os.path.realpath(__file__)):
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Anchor domains file:\t" + options.anchorfile)
    if options.anchorfile == "anchor_domains.txt":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("String for exclusion of gbk files:\t" + options.exclude_gbk_str)
    if options.exclude_gbk_str == "final":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("Skip hmmscan?:\t" + ("True" if options.skip_hmmscan else "False"))
    if not options.skip_hmmscan:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Skip MAFFT?:\t" + ("True" if options.skip_mafft else "False"))
    if not options.skip_mafft:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Skip all?:\t" + ("True" if options.skip_all else "False"))
    if not options.skip_all:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("Cutoff values for final network:\t" + options.sim_cutoffs)
    if options.sim_cutoffs == "1,0.85,0.75,0.6,0.4,0.2":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Neighborhood variable for GK:\t" + str(options.nbhood))
    if options.nbhood == 4:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
    
    pf.write("\nMAFFT parameters:\n")
    
    pf.write("Additional MAFFT parameters:\t\"" + options.mafft_pars + "\"")
    if options.mafft_pars == "":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    pf.write("Alignment method for MAFFT:\t\"" + options.al_method + "\"")
    if options.al_method == "--retree 2":
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    pf.write("Maxiterate:\t" + str(options.maxit))
    if options.maxit == 1000:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Threads in MAFFT:\t" + str(options.mafft_threads))
    if options.mafft_threads == -1:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")
        
    pf.write("Use own method for domain identity?:\t" + ("True" if options.use_perc_id else "False"))
    if options.use_perc_id:
        pf.write("\t(default)\n")
    else:
        pf.write("\n")

    pf.close()
