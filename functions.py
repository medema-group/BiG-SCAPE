#!/usr/bin/env python

"""
BiG-SCAPE

PI: Marnix Medema               marnix.medema@wur.nl

Maintainers:
Jorge Navarro                   jorge.navarromunoz@wur.nl

Developers:
Jorge Navarro                   jorge.navarromunoz@wur.nl
Satria Kautsar                  sakautsar@lbl.gov
Emmanuel (Emzo) de los Santos   E.De-Los-Santos@warwick.ac.uk

Functions used by bigscape.py


# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""

import os
import sys
import json

global verbose
verbose = False


def create_directory(path, kind, clean):
    # TODO consider makedirs(path,exists_ok=True)
    try:
        os.makedirs(path)
    except OSError as e:
        # 17 (Linux): "[Errno 17] File exists";
        # 183 (Windows) "[Error 183] Cannot create a file when that file already exists"
        if "Errno 17" in str(e) or "Error 183" in str(e):
            print(" " + kind + " folder already exists")
            if clean:
                print("  Cleaning folder")
                for thing in os.listdir(path):
                    os.remove(os.path.join(path,thing))
        else:
            print("Error: unexpected error creating " + kind + " folder")
            sys.exit(str(e))


def get_anchor_domains(filename):
    """Get the anchor/marker domains from a txt file.
    This text file should contain one Pfam id per line.
    A second column (separated by a tab) with more comments is allowed"""

    domains = set()
    
    try:
        with open(filename, "r") as handle:
            for line in handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    # ignore domain versions
                    domains.add(line.strip().split("\t")[0].split(".")[0])
        return domains
    except IOError:
        print("You have not provided the anchor_domains.txt file.")
        print("if you want to make use of the anchor domains in the DSS distance\
            metric, make a file that contains a Pfam domain on each line.")
        return set()
        

def get_domain_list(filename):
    """Convert the Pfam string in the .pfs files to a Pfam list"""
    
    domains = []
    with open(filename, "r") as handle:
        # note: cannot filter the domain version as the BGCs dictionary needs the
        # complete Pfam id
        domains = handle.readline().strip().split(" ")
        
    return domains


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
                    overlapping_aminoacids = overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4]))
                    overlap_perc_loc1 = overlap_perc(overlapping_aminoacids, int(row1[4])-int(row1[3]))
                    overlap_perc_loc2 = overlap_perc(overlapping_aminoacids, int(row2[4])-int(row2[3]))
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[1]) >= float(row2[1]): #see which has a better score
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
        
        # uses nucleotide coordinates for sorting
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
        # removing the domain version at this point is not a very good idea
        # because we need it complete for hmmalign
        domains.append(row[5]) #save the pfam domains for the .pfs file

    return pfd_matrix, domains                       
                        

def write_pfd(pfd_handle, matrix):
    for row in matrix:
        row = "\t".join(row)
        pfd_handle.write(row+"\n")
        
    pfd_handle.close() 
    

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

  
def BGC_dic_gen(filtered_matrix):
    """Generates the: { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } part of the BGCs variable."""
    
    # It would be tempting to do `bgc_dict = defaultdict(list)` but then the 
    # dictionary keeps the defaultdict definition. Later on in the distance
    # calculation function, we do a 
    #  ```try: unshared_occurrences = BGCs[A][unshared_domain]
    #  except KeyError: unshared_occurrences = BGCs[B][unshared_domain]```
    # Problem is the `try` would always be successful because of defaultdict and
    # an empty list is added for the unshared_domain!
    bgc_dict = {}
    for row in filtered_matrix:
        header = row[-1] + ":" + row[3] + ":" + row[4] # add domain positions
        try: #Should be faster than performing `if key in dictionary.keys()`
            bgc_dict[row[5]].append(header)
        except KeyError: # First occurrence of this domain in the cluster
           bgc_dict[row[5]] = [header]
            
    return bgc_dict

    
def save_domain_seqs(filtered_matrix, fasta_dict, domains_folder, outputbase):
    """Write fasta sequences for the domains in the right pfam-domain file"""
    for row in filtered_matrix:
        domain = row[5]
        header = row[-1].strip()
        seq = fasta_dict[header] #access the sequence by using the header
        
        domain_file = open(os.path.join(domains_folder, domain + ".fasta"), 'a') #append to existing file
        domain_file.write(">{}:{}:{}\n{}\n".format(header, row[3], row[4],
            seq[int(row[3])-1:int(row[4])])) #only use the range of the pfam domain within the sequence
        domain_file.close()


def write_network_matrix(matrix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info):
    """
    An entry in the distance matrix is currently (all floats):
      0         1       2      3      4    5    6    7    8        9    10    11        12          13      14
    clus1Idx clus2Idx  rawD  sqrtSim  Jac  DSS  AI rDSSna  rDSSa   S    Sa lcsStartA lcsStartB  seedLength reverse
    
    The final row in the network file is currently:
      0      1      2     3      4   5   6     7       8    9   10    11       12
    clus1  clus2  rawD  sqrtSim  J  DSS  AI  rDSSna  rDSSa  S   Sa  combGrp  ShrdGrp
    """
    
    #Open file handles for each cutoff
    networkfiles = {}
    cutoffs, filenames = zip(*cutoffs_and_filenames)
    for cutoff, filename in cutoffs_and_filenames:
        networkfiles[cutoff] = open(filename, "w")
        networkfiles[cutoff].write("Clustername 1\tClustername 2\tRaw distance\tSquared similarity\tJaccard index\tDSS index\tAdjacency index\traw DSS non-anchor\traw DSS anchor\tNon-anchor domains\tAnchor domains\tCombined group\tShared group\n")
      
    #Dictionaries to keep track of connected nodes, to know which are singletons
    clusterSetAllDict = {}
    clusterSetConnectedDict = {}
    for cutoff in cutoffs:
        clusterSetAllDict[cutoff] = set()
        clusterSetConnectedDict[cutoff] = set()

    for matrix_entry in matrix:
        gc1 = clusterNames[int(matrix_entry[0])]
        gc2 = clusterNames[int(matrix_entry[1])]
        row = [gc1, gc2]
        
        # get AntiSMASH annotations
        clus1group = bgc_info[gc1].product
        clus2group = bgc_info[gc2].product
        
        # add all the other floats
        row.extend(matrix_entry[2:-6])
        
        # add number of anchor/non-anchor domains as integers
        row.append(int(matrix_entry[-6]))
        row.append(int(matrix_entry[-5]))

        # prepare combined group
        if clus1group != "" and clus2group != "": #group1, group2
            row.append(" - ".join(sorted([clus1group,clus2group])))
        elif clus2group != "":
            row.append(clus2group)
        elif clus1group != "":
            row.append(clus1group)
        else:
            row.append("NA")
    
        # prepare share group (if they indeed share it)
        if clus1group == clus2group:
            row.append(clus1group)
        else:
            row.append("")

        for cutoff in cutoffs:
            clusterSetAllDict[cutoff].add(gc1)
            clusterSetAllDict[cutoff].add(gc2)
            
            if row[2] < cutoff:
                clusterSetConnectedDict[cutoff].add(gc1)
                clusterSetConnectedDict[cutoff].add(gc2)
                
                networkfiles[cutoff].write("\t".join(map(str,row)) + "\n")


    #Add the nodes without any edges, give them an edge to themselves with a distance of 0
    if include_singletons == True:
        for cutoff in cutoffs:
            for gc in clusterSetAllDict[cutoff]-clusterSetConnectedDict[cutoff]:
                #Arbitrary numbers for S and Sa domains: 1 of each (logical would be 0,0 but 
                # that could mess re-analysis with divisions-by-zero;
                networkfiles[cutoff].write("\t".join([gc, gc, "0", "1", "1", "1", "1", "0", "0", "1", "1", "", ""]) + "\n")

    #Close all files
    for networkfile in networkfiles.values():
        networkfile.close()


def fasta_parser(handle):
    """Parses a fasta file, and stores it in a dictionary.
    Only works if there are no duplicate fasta headers in the fasta file"""
    
    fasta_dict = {}
    header = ""
    for line in handle:
        if line[0] == ">":
            header=line.strip()[1:]
        else:
            try:
                fasta_dict[header] += line.strip()
            except KeyError:
                fasta_dict[header] = line.strip()

    return fasta_dict


def get_fasta_keys(handle):
    """Parses a fasta file, only stores headers
    """
    
    header_list = []
    for line in handle:
        if line[0] == ">":
            header_list.append(line.strip()[1:])
            
    return header_list


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


def sort_bgc(product):
    """Sort BGC by its type. Uses antiSMASH annotations
    (see 
    https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)
    """
    
    # TODO: according with current (2021-05) antiSMASH rules:
    # prodigiosin and PpyS-KS -> PKS
    # CDPS -> NRPS
    pks1_products = {'t1pks', 'T1PKS'}
    pksother_products = {'transatpks', 't2pks', 't3pks', 'otherks', 'hglks', 
                     'transAT-PKS', 'transAT-PKS-like', 'T2PKS', 'T3PKS', 
                     'PKS-like', 'hglE-KS', 'prodigiosin', 'PpyS-KS'}
    nrps_products = {'nrps', 'CDPS', 'NRPS', 'NRPS-like', 'thioamide-NRP', 
                    'NAPAA', 'mycosporine-like'}
    ripps_products = {'lantipeptide', 'thiopeptide', 'bacteriocin', 'linaridin', 
                     'cyanobactin', 'glycocin', 'LAP', 'lassopeptide', 
                     'sactipeptide', 'bottromycin', 'head_to_tail', 'microcin', 
                     'microviridin', 'proteusin', 'lanthipeptide', 'lipolanthine', 
                     'RaS-RiPP', 'fungal-RiPP', 'TfuA-related', 'guanidinotides', 
                     'RiPP-like', 'lanthipeptide-class-i', 'lanthipeptide-class-ii', 
                     'lanthipeptide-class-iii', 'lanthipeptide-class-iv',
                     'lanthipeptide-class-v', 'ranthipeptide', 'redox-cofactor',
                     'thioamitides', 'epipeptide', 'cyclic-lactone-autoinducer',
                     'spliceotide', 'RRE-containing', 'crocagin'}
    saccharide_products = {'amglyccycl', 'oligosaccharide', 'cf_saccharide', 
                     'saccharide'}
    others_products = {'acyl_amino_acids', 'arylpolyene', 'aminocoumarin', 
                       'ectoine', 'butyrolactone', 'nucleoside', 'melanin', 
                       'phosphoglycolipid', 'phenazine', 'phosphonate', 'other', 
                       'cf_putative', 'resorcinol', 'indole', 'ladderane', 
                       'PUFA', 'furan', 'hserlactone', 'fused', 'cf_fatty_acid', 
                       'siderophore', 'blactam', 'fatty_acid', 'PpyS-KS', 'CDPS', 
                       'betalactone', 'PBDE', 'tropodithietic-acid', 'NAGGN', 
                       'halogenated', 'pyrrolidine'}
    
    # PKS_Type I
    if product in pks1_products:
        return("PKSI")
    # PKS Other Types
    elif product in pksother_products:
        return("PKSother")
    # NRPs
    elif product in nrps_products:
        return("NRPS")
    # RiPPs
    elif product in ripps_products:
        return("RiPPs")
    # Saccharides
    elif product in saccharide_products:
        return("Saccharides")
    # Terpenes
    elif product == 'terpene':
        return("Terpene")
    # PKS/NRP hybrids
    elif len(product.split(".")) > 1:
        #print("  Possible hybrid: (" + cluster + "): " + product)
        # cf_fatty_acid category contains a trailing empty space
        subtypes = set(s.strip() for s in product.split("."))
        if len(subtypes - (pks1_products | pksother_products | nrps_products)) == 0:
            if len(subtypes - nrps_products) == 0:
                return("NRPS")
            elif len(subtypes - (pks1_products | pksother_products)) == 0:
                return("PKSother") # pks hybrids
            else:
                return("PKS-NRP_Hybrids")
        elif len(subtypes - ripps_products) == 0:
            return("RiPPs")
        elif len(subtypes - saccharide_products) == 0:
            return("Saccharide")
        else:
            return("Others") # other hybrid
    # Others
    elif product in others_products:
        return("Others")
    # ??
    elif product == "":
        # No product annotation. Perhaps not analyzed by antiSMASH
        return("Others")
    else:
        print("  Warning: unknown product '{}'".format(product))
        return("Others")


def write_parameters(output_folder, parameters):
    """Write a file with all the details of the run.
    Will overwrite previous versions"""
    
    with open(os.path.join(output_folder,"parameters.txt"), "w") as parameters_file:
        parameters_file.write(" ".join(parameters))
 


def generatePfamColorsMatrix(pfam_domain_colors):
    '''

    :param pfam_domain_colors: tab-delimited file
    :return: dictionary with pfam ID as key and rgb colors as value
    '''
    pfam_colors = {}

    if os.path.isfile(pfam_domain_colors):
        print("  Found file with Pfam domain colors")
        with open(pfam_domain_colors, "r") as cat_handle:
            for line in cat_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    domain = row[0]
                    rgb = row[-1]
                    pfam_colors[domain] = rgb
    else:
        print("  File pfam_domain_colors was NOT found")

    return pfam_colors


def add_to_bigscape_results_js(module_name, subs, result_js_file):
    bigscape_results = [];
    if os.path.isfile(result_js_file):
        with open(result_js_file, "r") as bs_js:
            line = bs_js.read()
            assert line.startswith("var bigscape_results = ")
            assert line.endswith(";")
            bigscape_results = json.loads(line[23:-1])
    bigscape_results.append({ "label" : module_name, "networks" : subs })
    with open(result_js_file, "w") as bs_js:
        bs_js.write("var bigscape_results = {};".format(json.dumps(bigscape_results, indent=4, separators=(',', ':'), sort_keys=True)))


def get_composite_bgc_similarities(bgcs_1, bgcs_2, sim_matrix):
    num_pairs = 0
    sum_sim = 0.00
    min_sim = (1.00, -1, -1)
    max_sim = (0.00, -1, -1)
    for bgc_1 in bgcs_1:
        for bgc_2 in bgcs_2:
            sim = 0.00 if (bgc_1 == bgc_2) else (sim_matrix[bgc_1][bgc_2] if ((bgc_1 in sim_matrix) and (bgc_2 in sim_matrix[bgc_1])) else sim_matrix[bgc_2][bgc_1])
            sum_sim += sim
            if (sim < min_sim[0]):
                min_sim = (sim, bgc_1, bgc_2) if (bgc_2 > bgc_1) else (sim, bgc_2, bgc_1)            
            if (sim > max_sim[0]):
                max_sim = (sim, bgc_1, bgc_2) if (bgc_2 > bgc_1) else (sim, bgc_2, bgc_1)
            num_pairs += 1
    return ((sum_sim / num_pairs), min_sim, max_sim) # let it throw infinite division error by itself
