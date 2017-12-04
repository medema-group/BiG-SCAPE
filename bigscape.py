#!/usr/bin/env python


"""
BiG-SCAPE

PI: Marnix Medema               marnix.medema@wur.nl

Developers:
Jorge Navarro                   jorge.navarromunoz@wur.nl
Emmanuel (Emzo) de los Santos   E.De-Los-Santos@warwick.ac.uk
Marley Yeong                    marleyyeong@live.nl


Dependencies: hmmer, biopython, (mafft), munkres.py

Usage:   Please see `python bigscape.py -h`

Example: python bigscape.py -c 8 --pfam_dir ./ -i ./inputfiles -o ./results

Status: development/testing

See more info on
https://git.wageningenur.nl/medema-group/BiG-SCAPE


# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""
# Makes sure the script can be used with Python 2 as well as Python 3.
from __future__ import print_function
from __future__ import division
from sys import version_info
if version_info[0]==2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0]==3:
    import pickle # for storing and retrieving dictionaries

from math import exp, log
import os
import subprocess
import sys
import time
from glob import glob
from itertools import combinations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from argparse import ArgumentParser
from difflib import SequenceMatcher
from operator import itemgetter

from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix
from Bio import Phylo

from functions import *
from ArrowerSVG import *

import numpy as np
from array import array
from scipy.sparse import lil_matrix
from scipy.optimize import linear_sum_assignment
import pysapc
import json
import shutil
from distutils import dir_util
from sklearn.cluster import AffinityPropagation

def get_gbk_files(inputdir, outputdir, bgc_fasta_folder, min_bgc_size, exclude_gbk_str, bgc_info):
    """Searches given directory for genbank files recursively, will assume that
    the genbank files that have the same name are the same genbank file. 
    Returns a dictionary that contains the names of the clusters found as keys
    and a list that contains [0] a path to the genbank file and [1] the 
    samples that the genbank file is a part of.
    Extract and write the sequences as fasta files if not already in the Fasta 
    folder.
    return: {cluster_name:[genbank_path,[s_a,s_b...]]}
    """

    genbankDict = {}

    file_counter = 0
    processed_sequences = 0
    biosynthetic_genes = set()
    product_list_per_record = []
    fasta_data = []
    save_fasta = True
    contig_edge = False
    total_seq_length = 0
    files_no_proteins = []
    files_no_biosynthetic_genes = []
    
    print("\nImporting GenBank files")
    if type(exclude_gbk_str) == str and exclude_gbk_str != "":
        print(" Skipping files with '{}' in their filename".format(exclude_gbk_str))
    elif type(exclude_gbk_str) == list and exclude_gbk_str != []:
        print(" Skipping files with one or more of the following strings in \
            their filename: {}".format(", ".join(exclude_gbk_str)))

    current_dir = ""
    for dirpath, dirnames, filenames in os.walk(inputdir):
        head, tail = os.path.split(dirpath)

        if current_dir != tail:
            current_dir = tail

        genbankfilelist = []

        for fname in filenames:
            if fname[-3:] != "gbk":
                continue
            
            clusterName = fname[:-4]
            
            if type(exclude_gbk_str) == str and exclude_gbk_str != "" and \
                                                exclude_gbk_str in fname:
                print(" Skipping file " + fname)
                continue
            elif type(exclude_gbk_str) == list and exclude_gbk_str != [] and \
                            any([word in fname for word in exclude_gbk_str]):
                print(" Skipping file " + fname)
                continue
            if "_ORF" in fname:
                print(" Skipping file {} (string '_ORF' is used internally)".format(fname))
                continue
            
            if " " in fname:
                sys.exit("\nError: Input GenBank files should not have spaces in their filenames as HMMscan cannot process them properly ('too many arguments').")
            
            # See if we need to keep the sequence
            outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                save_fasta = False
            else:
                save_fasta = True
                
            try:
                # basic file verification. Substitutes check_data_integrity
                records = list(SeqIO.parse(os.path.join(dirpath,fname), "genbank"))
            except ValueError as e:
                print("   Error with file {}: \n    '{}'".format(os.path.join(dirpath, fname), str(e)))
                print("    (This file will be excluded from the analysis)")
                continue
            else:
                total_seq_length = 0
                bgc_size = 0
                cds_ctr = 0
                product = "no type"
                del product_list_per_record[:]
                
                max_width = 0 # This will be used for the SVG figure
                record_count = 0
                
                for record in records:
                    record_count += 1
                    bgc_size += len(record.seq)
                    if len(record.seq) > max_width:
                        max_width = len(record.seq)
                    
                    for feature in record.features:
                        if "cluster" in feature.type:
                            if "product" in feature.qualifiers:
                                if len(feature.qualifiers["product"]) > 1:
                                    print("  WARNING: more than product annotated in record " + str(record_count) + ", " + fname)
                                    break
                                else:
                                    product_list_per_record.append(feature.qualifiers["product"][0].replace(" ",""))
                            if "contig_edge" in feature.qualifiers:
                                # there might be mixed contig_edge annotations
                                # in multi-record files. Turn on contig_edge when
                                # there's at least one annotation
                                if feature.qualifiers["contig_edge"][0] == "True":
                                    if verbose:
                                        print(" Contig edge detected in {}".format(fname))
                                    contig_edge = True
                
                        # Get biosynthetic genes + sequences
                        if feature.type == "CDS":
                            cds_ctr += 1
                                    
                            CDS = feature
                            gene_id = ""
                            if "gene" in CDS.qualifiers:
                                # In principle, we should keep a list of genes with
                                # the same id (isoforms) and only keep the largest
                                # TODO
                                gene_id = CDS.qualifiers.get('gene',"")[0]
                            
                            protein_id = ""
                            if "protein_id" in CDS.qualifiers:
                                protein_id = CDS.qualifiers.get('protein_id',"")[0]
                            
                            # nofuzzy_start/nofuzzy_end are obsolete
                            # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#nofuzzy_start
                            gene_start = max(0, int(CDS.location.start))
                            gene_end = max(0, int(CDS.location.end))
                            direction = CDS.location.strand
                            
                            if direction == 1:
                                strand = '+'
                            else:
                                strand = '-'
                            fasta_header = "{}_ORF{}:gid:{}:pid:{}:loc:{}:{}:strand:{}".format(clusterName, str(cds_ctr), str(gene_id), str(protein_id), str(gene_start), str(gene_end), strand)
                            fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this
                            fasta_header = fasta_header.replace(" ", "") #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header

                            if "sec_met" in feature.qualifiers:
                                if "Kind: biosynthetic" in feature.qualifiers["sec_met"]:
                                    biosynthetic_genes.add(fasta_header)

                            fasta_header = ">"+fasta_header
                            

                            if 'translation' in CDS.qualifiers.keys():
                                prot_seq = CDS.qualifiers['translation'][0]
                            # If translation isn't available translate manually, this will take longer
                            else:
                                nt_seq = CDS.location.extract(record.seq)
                                
                                # If we know sequence is an ORF (like all CDSs), codon table can be
                                #  used to correctly translate alternative start codons.
                                #  see http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25
                                # If the sequence has a fuzzy start/end, it might not be complete,
                                # (therefore it might not be the true start codon)
                                # However, in this case, if 'translation' not available, assume 
                                #  this is just a random sequence 
                                complete_cds = False 
                                
                                # More about fuzzy positions
                                # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc39
                                fuzzy_start = False 
                                if str(CDS.location.start)[0] in "<>":
                                    complete_cds = False
                                    fuzzy_start = True
                                    
                                fuzzy_end = False
                                if str(CDS.location.end)[0] in "<>":
                                    fuzzy_end = True
                                
                                #for protein sequence if it is at the start of the entry assume 
                                # that end of sequence is in frame and trim from the beginning
                                #if it is at the end of the genbank entry assume that the start 
                                # of the sequence is in frame
                                reminder = len(nt_seq)%3
                                if reminder > 0:
                                    if fuzzy_start and fuzzy_end:
                                        print("Warning, CDS ({}, {}) has fuzzy start and end positions, and a sequence length not multiple of three. Skipping".format(clusterName, CDS.qualifiers.get('locus_tag',"")[0]))
                                        break
                                    
                                    if fuzzy_start:
                                        if reminder == 1:
                                            nt_seq = nt_seq[1:]
                                        else:
                                            nt_seq = nt_seq[2:]
                                    # fuzzy end
                                    else:
                                        #same logic reverse direction
                                        if reminder == 1:
                                            nt_seq = nt_seq[:-1]
                                        else:
                                            nt_seq = nt_seq[:-2]
                                
                                # The Genetic Codes: www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
                                if "transl_table" in CDS.qualifiers.keys():
                                    CDStable = CDS.qualifiers.get("transl_table", "")[0]
                                    prot_seq = str(nt_seq.translate(table=CDStable, to_stop=True, cds=complete_cds))
                                else:
                                    prot_seq = str(nt_seq.translate(to_stop=True, cds=complete_cds))
                                    
                            total_seq_length += len(prot_seq)
                            if save_fasta:
                                fasta_data.append((fasta_header, prot_seq))
        
                    # TODO: if len(biosynthetic_genes) == 0, traverse record again
                    # and add CDS with genes that contain domains labeled sec_met
                    # we'll probably have to have a list of domains if we allow
                    # fasta files as input
        
                if bgc_size > min_bgc_size:  # exclude the bgc if it's too small
                    file_counter += 1
                    # check what we have product-wise
                    # In particular, handle different products for multi-record files
                    product_set = set(product_list_per_record)
                    if len(product_set) == 1: # only one type of product
                        product = product_list_per_record[0]
                    elif "other" in product_set: # more than one, and it contains "other"
                        if len(product_set) == 2:
                            product = list(product_set - set(['other']))[0] # product = not "other"
                        else:
                            product = "-".join(product_set - set(['other'])) # likely a hybrid
                    else:
                        product = "-".join(product_set) # likely a hybrid
                     
                    
                    # assuming that the definition field is the same in all records
                    # product: antiSMASH predicted class of metabolite
                    # gbk definition
                    # number of records (for Arrower figures)
                    # max_width: width of the largest record (for Arrower figures)
                    # id: the GenBank's accession
                    #bgc_info[clusterName] = (product, records[0].description, len(records), max_width, records[0].id, biosynthetic_genes.copy())
                    # TODO contig_edge annotation is not present for antiSMASH v < 4
                    # Perhaps we can try to infer if it's in a contig edge: if
                    # - first biosynthetic gene start < 10kb or
                    # - max_width - last biosynthetic gene end < 10kb (but this will work only for the largest record)
                    bgc_info[clusterName] = bgc_data(records[0].id, records[0].description, product, len(records), max_width, records[0].annotations["organism"], ",".join(records[0].annotations["taxonomy"]), biosynthetic_genes.copy(), contig_edge)

                    if len(bgc_info[clusterName].biosynthetic_genes) == 0:
                        files_no_biosynthetic_genes.append(fname)

                    # TODO why re-process everything if it was already in the list?
                    # if name already in genbankDict.keys -> add current_dir
                    # else: extract all info
                    if clusterName in genbankDict.keys():
                        # Name was already in use. Use current_dir as the new sample's name
                        genbankDict[clusterName][1].add(current_dir) 
                    else:
                        # See if we need to write down the sequence
                        if total_seq_length > 0:
                            # location of first instance of the file is genbankDict[clustername][0]
                            genbankDict.setdefault(clusterName, [os.path.join(dirpath, fname), set([current_dir])])

                            if save_fasta:
                                processed_sequences += 1
                                with open(outputfile,'w') as fastaHandle:
                                    for header_sequence in fasta_data:
                                        fastaHandle.write("{}\n".format(str(header_sequence[0])))
                                        fastaHandle.write("{}\n".format(str(header_sequence[1])))
                        else:
                            files_no_proteins.append(fname)

                    if verbose:
                        print("  Adding {} ({} bps)".format(fname, str(bgc_size)))
                        
                else:
                    print(" Discarding {} (size less than {} bp, was {})".format(clusterName, str(min_bgc_size), str(bgc_size)))
                
                del fasta_data[:]
                biosynthetic_genes.clear()
    
    if file_counter == 0:
        sys.exit("\nError: There are no files to process")
        
    if file_counter == 1:
        sys.exit("\nError: Only one file found. Please input at least two files")
        
    if len(files_no_proteins) > 0:
        print("  Warning: Input set has files without protein sequences. They will be discarded")
        print("   (See no_sequences_list.txt)")
        with open(os.path.join(outputdir, "no_sequences_list.txt"), "w") as noseqs:
            for f in sorted(files_no_proteins):
                noseqs.write("{}\n".format(f))
        
    if len(files_no_biosynthetic_genes) > 0 and (mode == "lcs" or mode == "local"):
        print("  Warning: Input set has files with no Biosynthetic Genes (affects alignment mode)")
        print("   See no_biosynthetic_genes_list.txt")
        with open(os.path.join(outputdir, "logs", "no_biosynthetic_genes_list.txt"), "w") as nobiogenes:
            for f in sorted(files_no_biosynthetic_genes):
                nobiogenes.write("{}\n".format(f))
    
    print("\n Starting with {:d} files".format(file_counter))
    print(" Files that had its sequence extracted: {:d}".format(processed_sequences))

    return genbankDict


def timeit(funct):
    """Writes the runtimes of functions to a file and prints them on screen.

    Input:
    - funct: a function
    Output:
    - That function its output.
    - Runtime of that function in a file called commands.txt and on screen.
    """
    def _wrap(*args):
        start_time = time.time()
        ret = funct(*args)
        runtime = time.time()-start_time
        runtime_string = '{} took {:.3f} seconds'.format(funct.__name__, runtime)
        # To prevent insignificant runtimes from ending up in the file.
        if runtime > 1:
            with open(os.path.join(log_folder, "runtimes.txt"), 'a') as timings_file:
                timings_file.write(runtime_string + "\n")
            print(runtime_string)
        return ret
    
    return _wrap


@timeit
def generate_network(cluster_pairs, cores):
    """Distributes the distance calculation part
    cluster_pairs is a list of triads (cluster1_index, cluster2_index, BGC class)
    """
    
    pool = Pool(cores, maxtasksperchild=100)
    
    #Assigns the data to the different workers and pools the results back into
    # the network_matrix variable
    network_matrix = pool.map(generate_dist_matrix, cluster_pairs)

    # --- Serialized version of distance calculation ---
    # For the time being, use this if you have memory issues
    #network_matrix = []
    #for pair in cluster_pairs:
      #network_matrix.append(generate_dist_matrix(pair))

    return network_matrix


def generate_dist_matrix(parms):
    """Unpack data to actually launch cluster_distance for one pair of BGCs"""
    
    cluster1Idx,cluster2Idx,bgcClassIdx = [int(parm) for parm in parms]
    cluster1 = clusterNames[cluster1Idx]
    cluster2 = clusterNames[cluster2Idx]
    bgc_class = bgcClassNames[bgcClassIdx]

    try:
        domain_list_A = DomainList[cluster1]
        domain_list_B = DomainList[cluster2]
    except KeyError:
        print(" Warning: domain list for {} or {} was not found. Extracting from pfs files".format(cluster1, cluster2))
        
        cluster_file1 = os.path.join(output_folder, cluster1 + ".pfs")
        cluster_file2 = os.path.join(output_folder, cluster2 + ".pfs")
        
        domain_list_A = get_domain_list(cluster_file1)
        domain_list_B = get_domain_list(cluster_file2)
    
    # this really shouldn't happen if we've filtered domain-less gene clusters already
    if len(domain_list_A) == 0 or len(domain_list_B) == 0:
        print("   Warning: Regarding distance between clusters {} and {}:".format(cluster1, cluster2))
        if len(domain_list_A) == 0 and len(domain_list_B) == 0:
            print("   None have identified domains. Distance cannot be calculated")
        elif (domain_list_A) == 0:            
            print("   Cluster {} has no identified domains. Distance set to 1".format(cluster1))
        else:
            print("   Cluster {} has no identified domains. Distance set to 1".format(cluster2))
            
        # last two values (S, Sa) should really be zero but this could give rise to errors when parsing 
        # the network file (unless we catched the case S = Sa = 0

        # cluster1Idx, cluster2Idx, distance, jaccard, DSS, AI, rDSSNa, rDSSa, 
        #   S, Sa, lcsStartA, lcsStartB
        return array('f',[cluster1Idx,cluster2Idx,1,0,0,0,0,0,1,1,0,0])
    
    # "Domain Count per Gene". List of simple labels (integers) indicating number
    # of domains belonging to each gene
    dcg_a = DomainCountGene[cluster1]
    dcg_b = DomainCountGene[cluster2]
    
    # Position of the anchor genes (i.e. genes with domains in the anchor
    # domain list). Should probably be the Core Biosynthetic genes marked by
    # antiSMASH
    core_pos_a = corebiosynthetic_position[cluster1]
    core_pos_b = corebiosynthetic_position[cluster2]
    
    # go = "gene orientation"
    go_a = BGCGeneOrientation[cluster1]
    go_b = BGCGeneOrientation[cluster2]
    
    dist, jaccard, dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB, seedLength, reverse = cluster_distance_lcs(cluster1, 
            cluster2, domain_list_A, domain_list_B, dcg_a, dcg_b, core_pos_a, 
            core_pos_b, go_a, go_b, bgc_class)
        
    network_row = array('f',[cluster1Idx, cluster2Idx, dist, (1-dist)**2, jaccard, 
                             dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB, seedLength, reverse])
    return network_row
    

def score_expansion(x_string_, y_string_, downstream):
    """
    Input:
    x_string: list of strings. This one tries to expand based on max score
    y_string: reference list. This was already expanded.
    downstream: If true, expansion goes from left to right.
    
    Output:
    max_score, a
    Where a is length of expansion.
    """
    match = 5
    mismatch = -3
    gap = -2
        
    if downstream:
        x_string = x_string_
        y_string = y_string_
    else:
        # Expansion goes upstream. It's easier to flip both slices and proceed
        # as if going downstream.
        x_string = list(reversed(x_string_))
        y_string = list(reversed(y_string_))
    

    # how many gaps to open before calling a mismatch is more convenient?
    max_gaps_before_mismatch = abs(int(match/gap))
    max_score = 0
    score = 0

    pos_y = 0
    a = 0
    b = 0
    for pos_x in range(len(x_string)):
        g = x_string[pos_x]
        
        # try to find g within the rest of the slice
        # This has the obvious problem of what to do if a gene is found _before_
        # the current y_slice (translocation). Duplications could also complicate
        # things
        try:
            match_pos = y_string[pos_y:].index(g)
        except ValueError:
            score += mismatch
        else:
            score += match + match_pos*gap
            pos_y += match_pos + 1 # move pointer one position after match
            
            # Greater or equals. 'equals' because even if max_score isn't 
            # larger, it's an expansion
            if score >= max_score:
                max_score = score
                # keep track of expansion. Account for zero-based numbering
                a = pos_x + 1
                            
        #print(pos_x, score, status, g)
        #print("")
            
    return max_score, a


def cluster_distance_lcs(A, B, A_domlist, B_domlist, dcg_A, dcg_b, core_pos_A, core_pos_b, go_A, go_b, bgc_class):
    """Compare two clusters using information on their domains, and the 
    sequences of the domains. 
    This version first tries to search for the largest common slices of both BGCs by 
    finding the Longest Common Substring of genes (based on domain content), also
    trying the reverse on one of the BGCs (capital "B" will be used when the 'true'
    orientation of the BGC is found)
    Then the algorithm will try to expand the slices to either side. Once each 
    slice are found, they have to pass one more validation: a core biosynthetic 
    gene (marked in antiSMASH) must be present.
    Capital letters indicate the "true" orientation of the BGC (first BGC, 'A'
    is kept fixed)
    
    Output:
    raw distance, jaccard index, domain sequence similarity, adjacency index, 
    raw DSS non-anchor, raw DSS anchor, number of non-anchor domains in the pair
    (S), number of anchor domains in the pair
    
    """
    
    Jaccardw, DSSw, AIw, anchorboost = bgc_class_weight[bgc_class]

    temp_domain_fastas = {}
    
    # Number of genes in each BGC
    lenG_A = len(dcg_A)
    lenG_B = len(dcg_b)
    
    setA = set(A_domlist)
    setB = set(B_domlist)
    intersect = setA & setB
    
    S = 0
    S_anchor = 0
    
    # define the subset of domain sequence tags to include in
    # the DSS calculation. This is done for every domain.
    A_domain_sequence_slice_bottom = defaultdict(int)
    A_domain_sequence_slice_top = defaultdict(int)
    B_domain_sequence_slice_bottom = defaultdict(int)
    B_domain_sequence_slice_top = defaultdict(int)
    
    # Detect totally unrelated pairs from the beginning
    if len(intersect) == 0:
        # Count total number of anchor and non-anchor domain to report in the 
        # network file. Apart from that, these BGCs are totally unrelated.
        for domain in setA:
            # This is a bit of a hack. If pfam domain ids ever change in size
            # we'd be in trouble. The previous approach was to .split(".")[0]
            # but it's more costly
            if domain[:7] in anchor_domains:
                S_anchor += len(BGCs[A][domain])
            else:
                S += len(BGCs[A][domain])
                
        for domain in setB:
            if domain[:7] in anchor_domains:
                S_anchor += len(BGCs[B][domain])
            else:
                S += len(BGCs[B][domain])
        
        return 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, S, S_anchor, 0, 0


    # initialize domain sequence slices
    # They might change if we manage to find a valid overlap
    for domain in setA:
        A_domain_sequence_slice_bottom[domain] = 0
        A_domain_sequence_slice_top[domain] = len(BGCs[A][domain])
        
    for domain in setB:
        B_domain_sequence_slice_bottom[domain] = 0
        B_domain_sequence_slice_top[domain] = len(BGCs[B][domain])
        
    #initialize domlist borders for AI
    domA_start = 0
    domA_end = len(A_domlist)
    domB_start = 0
    domB_end = len(B_domlist)

    # always find lcs seed to use for offset alignment in visualization
    
    # Compress the list of domains according to gene information. For example:
    # A_domlist = a b c d e f g
    # dcg_a =   1  3  1  2 Number of domains per each gene in the BGC
    # go_a =    1 -1 -1  1 Orientation of each gene
    # A_string = a dcb e fg List of concatenated domains
    # Takes into account gene orientation. This works effectively as putting all
    # genes in the same direction in order to be able to compare their domain content
    A_string = []
    start = 0
    for g in range(lenG_A):
        domain_count = dcg_A[g]
        if go_A[g] == 1:
            # x[2:] <- small optimization, drop the "PF" from the pfam ids
            A_string.append("".join(x[2:] for x in A_domlist[start:start+domain_count]))
        else:
            A_string.append("".join(A_domlist[x][2:] for x in range(start+domain_count-1, start-1 ,-1)))
        start += domain_count
        
    b_string = []
    start = 0
    for g in range(lenG_B):
        domain_count = dcg_b[g]
        if go_b[g] == 1:
            b_string.append("".join(x[2:] for x in B_domlist[start:start+domain_count]))
        else:
            b_string.append("".join(B_domlist[x][2:] for x in range(start+domain_count-1, start-1 ,-1)))
        start += domain_count
        
    b_string_reverse = list(reversed(b_string))
        
    seqmatch = SequenceMatcher(None, A_string, b_string)
    # a: start position in A
    # b: start position in B
    # s: length of the match
    a, b, s = seqmatch.find_longest_match(0, lenG_A, 0, lenG_B)
    #print(A, B)
    #print(a, b, s)
    
    seqmatch = SequenceMatcher(None, A_string, b_string_reverse)
    ar, br, sr = seqmatch.find_longest_match(0, lenG_A, 0, lenG_B)
    #print(ar, br, sr)
    
    # We need to keep working with the correct orientation
    if s >= sr:
        dcg_B = dcg_b
        B_string = b_string
        # note: these slices are in terms of genes, not domains (which are 
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        sliceStartA = a
        sliceStartB = b
        sliceLengthA = s
        sliceLengthB = s
        
        reverse = False
        b_name = B
        
    else:
        dcg_B = list(reversed(dcg_b))
        B_string = b_string_reverse
        
        sliceStartA = ar
        sliceStartB = br
        sliceLengthA = sr
        sliceLengthB = sr
        
        # We'll need to know if we're working in reverse so that the start 
        # postion of the final slice can be transformed to the original orientation
        reverse = True
        b_name = B + "*"
        
    lcsStartA = sliceStartA
    lcsStartB = sliceStartB
    seedLength = sliceLengthA
    
    
    # paint stuff on screen
    #if sliceStartB > sliceStartA:
        #offset_A = sliceStartB - sliceStartA
        #offset_B = 0
    #else:
        #offset_A = 0
        #offset_B = sliceStartA - sliceStartB
    #print(A, B)
    #print("  "*offset_A + " ".join(map(str,dcg_A[:sliceStartA])) + "[" + " ".join(map(str,dcg_A[sliceStartA:sliceStartA+sliceLengthA])) + "]" + " ".join(map(str,dcg_A[sliceStartA+sliceLengthA:])) + "\t" + A )
    #print("  "*offset_B + " ".join(map(str,dcg_B[:sliceStartB])) + "[" + " ".join(map(str,dcg_B[sliceStartB:sliceStartB+sliceLengthB])) + "]" + " ".join(map(str,dcg_B[sliceStartB+sliceLengthB:])) + "\t" + b_name)
    ##print(sliceStartA, sliceStartB, sliceLengthA)
    #print("")

    if mode=="lcs" or (mode=="auto" and (bgc_info[A].contig_edge or bgc_info[B].contig_edge)):
        #X: bgc that drive expansion
        #Y: the other bgc
        # forward: True if expansion is to the right
        # returns max_score, final positions for X and Y
        #score_expansion(X_string, dcg_X, Y_string, dcg_Y, downstream=True/False)
            
        # Expansion is relatively costly. We ask for a minimum of 3 genes
        # for the core overlap before proceeding with expansion.
        if sliceLengthA >= 3:
            # LEFT SIDE
            # Find which bgc has the least genes to the left. If both have the same 
            # number, find the one that drives the expansion with highest possible score
            if sliceStartA == sliceStartB:
                # assume complete expansion of A, try to expand B
                score_B, sbB = score_expansion(B_string[:sliceStartB], A_string[:sliceStartA], False)
                # assume complete expansion of B, try to expand A
                score_A, saA = score_expansion(A_string[:sliceStartA], B_string[:sliceStartB], False)
                
                if score_A > score_B or (score_A == score_B and saA > sbB):
                    sliceLengthA += saA
                    sliceLengthB += len(B_string[:sliceStartB])
                    sliceStartA -= saA
                    sliceStartB = 0
                else:
                    sliceLengthA += len(A_string[:sliceStartA])
                    sliceLengthB += sbB
                    sliceStartA = 0
                    sliceStartB -= sbB

            else:
                # A is shorter upstream. Assume complete extension. Find B's extension
                if sliceStartA < sliceStartB:
                    score_B, sb = score_expansion(B_string[:sliceStartB], A_string[:sliceStartA], False)
                    
                    sliceLengthA += len(A_string[:sliceStartA])
                    sliceLengthB += sb
                    sliceStartA = 0
                    sliceStartB -= sb
                else:
                    score_A, sa = score_expansion(A_string[:sliceStartA], B_string[:sliceStartB], False)
                    
                    sliceLengthA += sa
                    sliceLengthB += len(B_string[:sliceStartB])
                    sliceStartA -= sa
                    sliceStartB = 0

            # RIGHT SIDE
            # check which side is the shortest downstream. If both BGCs have the same
            # length left, choose the one with the best expansion
            downstream_A = lenG_A - sliceStartA - sliceLengthA
            downstream_B = lenG_B - sliceStartB - sliceLengthB
            if downstream_A == downstream_B:
                # assume complete extension of A, try to expand B
                score_B, xb = score_expansion(B_string[sliceStartB+sliceLengthB:], A_string[sliceStartA+sliceLengthA:], True)
                # assume complete extension of B, try to expand A
                score_A, xa = score_expansion(A_string[sliceStartA+sliceLengthA:], B_string[sliceStartB+sliceLengthB:], True)
                
            
                if (score_A == score_B and xa > xb) or score_A > score_B:
                    sliceLengthA += xa
                    sliceLengthB += len(B_string[sliceStartB+sliceLengthB:])
                else:
                    sliceLengthA += len(A_string[sliceStartA+sliceLengthA:])
                    sliceLengthB += xb
                    
            else:
                if downstream_A < downstream_B:
                    # extend all of remaining A
                    score_B, xb = score_expansion(B_string[sliceStartB+sliceLengthB:], A_string[sliceStartA+sliceLengthA:], True)
                    
                    sliceLengthA += len(A_string[sliceStartA+sliceLengthA:])
                    sliceLengthB += xb              
                    
                else:
                    score_A, xa = score_expansion(A_string[sliceStartA+sliceLengthA:], B_string[sliceStartB+sliceLengthB:], True)
                    sliceLengthA += xa
                    sliceLengthB += len(B_string[sliceStartB+sliceLengthB:])
        
        #print("  "*offset_A + " ".join(map(str,dcg_A[:sliceStartA])) + "[" + " ".join(map(str,dcg_A[sliceStartA:sliceStartA+sliceLengthA])) + "]" + " ".join(map(str,dcg_A[sliceStartA+sliceLengthA:])) + "\t" + A )
        #print("  "*offset_B + " ".join(map(str,dcg_B[:sliceStartB])) + "[" + " ".join(map(str,dcg_B[sliceStartB:sliceStartB+sliceLengthB])) + "]" + " ".join(map(str,dcg_B[sliceStartB+sliceLengthB:])) + "\t" + b_name)

        #print()
        #print(A_string)
        #print(B_string)
        #print()
        if min(sliceLengthA, sliceLengthB) >= 5:
            # First test passed. Find if there is a biosynthetic gene in both slices
            # (note that even if they are, currently we don't check whether it's 
            # actually the _same_ gene)
            biosynthetic_hit_A = False
            biosynthetic_hit_B = False
            
            for biosynthetic_position in core_pos_A:
                if biosynthetic_position >= sliceStartA and biosynthetic_position <= (sliceStartA+sliceLengthA):
                    biosynthetic_hit_A = True
                    break
            
            # return to original orientation if needed
            if reverse:
                sliceStartB = lenG_B - sliceStartB - sliceLengthB
                
            # using original core_pos_b
            for biosynthetic_position in core_pos_b:
                if biosynthetic_position >= sliceStartB and biosynthetic_position <= (sliceStartB + sliceLengthB):
                    biosynthetic_hit_B = True
                    break
                
            # finally...
            if biosynthetic_hit_A and biosynthetic_hit_B:
                domA_start = sum(dcg_A[:sliceStartA])
                domA_end = domA_start + sum(dcg_A[sliceStartA:sliceStartA+sliceLengthA])
                setA = set(A_domlist[domA_start:domA_end])
                
                domB_start = sum(dcg_b[:sliceStartB])
                domB_end = domB_start + sum(dcg_b[sliceStartB:sliceStartB+sliceLengthB])
                setB = set(B_domlist[domB_start:domB_end])
                
                intersect = setA & setB
                
                # re-adjust the indices for each domain so we get only the sequence
                # tags in the selected slice. First step: find out which is the 
                # first copy of each domain we're using
                for domain in A_domlist[:domA_start]:
                    A_domain_sequence_slice_bottom[domain] += 1
                for domain in B_domlist[:domB_start]:
                    B_domain_sequence_slice_bottom[domain] += 1
                    
                # Step 2: work with the last copy of each domain. 
                # Step 2a: make top = bottom
                for domain in setA:
                    A_domain_sequence_slice_top[domain] = A_domain_sequence_slice_bottom[domain]
                for domain in setB:
                    B_domain_sequence_slice_top[domain] = B_domain_sequence_slice_bottom[domain]
                
                # Step 2b: increase top with the domains in the slice
                for domain in A_domlist[domA_start:domA_end]:
                    A_domain_sequence_slice_top[domain] += 1
                for domain in B_domlist[domB_start:domB_end]:
                    B_domain_sequence_slice_top[domain] += 1
                    
            #else:
                #print(" - - Not a valid overlap found - - (no biosynthetic genes)\n")
                
        #else:
                #print(" - - Not a valid overlap found - - (shortest slice not large enough)\n")
                    
    # JACCARD INDEX
    Jaccard = len(intersect) / len(setA | setB)


    # DSS INDEX
    #domain_difference: Difference in sequence per domain. If one cluster does
    # not have a domain at all, but the other does, this is a (complete) 
    # difference in sequence 1. If both clusters contain the domain once, and 
    # the sequence is the same, there is a seq diff of 0.
    #S: Max occurence of each domain
    domain_difference_anchor,S_anchor = 0,0
    domain_difference,S = 0,0
        
    not_intersect = setA.symmetric_difference(setB)
        
    # Case 1
    #no need to look at seq identity, since these domains are unshared
    for unshared_domain in not_intersect:
        #for each occurence of an unshared domain do domain_difference += count 
        # of domain and S += count of domain
        unshared_occurrences = []

        try:
            num_unshared = A_domain_sequence_slice_top[unshared_domain] - A_domain_sequence_slice_bottom[unshared_domain]
        except KeyError:
            num_unshared = B_domain_sequence_slice_top[unshared_domain] - B_domain_sequence_slice_bottom[unshared_domain]
            
        # don't look at domain version, hence the split
        if unshared_domain[:7] in anchor_domains:
            domain_difference_anchor += num_unshared
        else:
            domain_difference += num_unshared
                    
    S = domain_difference # can be done because it's the first use of these
    S_anchor = domain_difference_anchor
        
    # Cases 2 and 3 (now merged)
    missing_aligned_domain_files = []
    for shared_domain in intersect:
        specific_domain_list_A = BGCs[A][shared_domain]
        specific_domain_list_B = BGCs[B][shared_domain]
        
        num_copies_a = A_domain_sequence_slice_top[shared_domain] - A_domain_sequence_slice_bottom[shared_domain]
        num_copies_b = B_domain_sequence_slice_top[shared_domain] - B_domain_sequence_slice_bottom[shared_domain]
        
        temp_domain_fastas.clear()
        
        accumulated_distance = 0
            
        # Fill distance matrix between domain's A and B versions
        DistanceMatrix = np.ndarray((num_copies_a,num_copies_b))
        for domsa in range(num_copies_a):
            for domsb in range(num_copies_b):
                sequence_tag_a = specific_domain_list_A[domsa + A_domain_sequence_slice_bottom[shared_domain]]
                sequence_tag_b = specific_domain_list_B[domsb + B_domain_sequence_slice_bottom[shared_domain]]
                
                seq_length = 0
                matches = 0
                gaps = 0
                
                try:
                    aligned_seqA = AlignedDomainSequences[sequence_tag_a]
                    aligned_seqB = AlignedDomainSequences[sequence_tag_b]
                    
                except KeyError:
                    # For some reason we don't have the multiple alignment files. 
                    # Try manual alignment
                    if shared_domain not in missing_aligned_domain_files and verbose:
                        # this will print everytime an unfound <domain>.algn is not found for every
                        # distance calculation (but at least, not for every domain pair!)
                        print("  Warning: {}.algn not found. Trying pairwise alignment...".format(shared_domain))
                        missing_aligned_domain_files.append(shared_domain)
                    
                    try:
                        unaligned_seqA = temp_domain_fastas[sequence_tag_a]
                        unaligned_seqB = temp_domain_fastas[sequence_tag_b]
                    except KeyError:
                        # parse the file for the first time and load all the sequences
                        with open(os.path.join(domains_folder, shared_domain + ".fasta"),"r") as domain_fasta_handle:
                            temp_domain_fastas = fasta_parser(domain_fasta_handle)
                        
                        unaligned_seqA = temp_domain_fastas[sequence_tag_a]
                        unaligned_seqB = temp_domain_fastas[sequence_tag_b]
                        
                    # gap_open = -15
                    # gap_extend = -6.67. These parameters were set up by Emzo
                    alignScore = pairwise2.align.globalds(unaligned_seqA, unaligned_seqB, scoring_matrix, -15, -6.67, one_alignment_only=True)
                    bestAlignment = alignScore[0]
                    aligned_seqA = bestAlignment[0]
                    aligned_seqB = bestAlignment[1]
                    
                    
                # - Calculate aligned domain sequences similarity -
                # Sequences *should* be of the same length unless something went
                # wrong elsewhere
                if len(aligned_seqA) != len(aligned_seqB):
                    print("\tWARNING: mismatch in sequences' lengths while calculating sequence identity ({})".format(shared_domain))
                    print("\t  Specific domain 1: {} len: {}".format(sequence_tag_a, str(len(aligned_seqA))))
                    print("\t  Specific domain 2: {} len: {}".format(sequence_tag_b, str(len(aligned_seqB))))
                    seq_length = min(len(aligned_seqA), len(aligned_seqB))
                else:
                    seq_length = len(aligned_seqA)
                    
                for position in range(seq_length):
                    if aligned_seqA[position] == aligned_seqB[position]:
                        if aligned_seqA[position] != "-":
                            matches += 1
                        else:
                            gaps += 1
                            
                DistanceMatrix[domsa][domsb] = 1 - ( matches/(seq_length-gaps) )
                
        #Only use the best scoring pairs
        BestIndexes = linear_sum_assignment(DistanceMatrix)
        accumulated_distance = DistanceMatrix[BestIndexes].sum()
        
        # the difference in number of domains accounts for the "lost" (or not duplicated) domains
        sum_seq_dist = (abs(num_copies_a-num_copies_b) + accumulated_distance)  #essentially 1-sim
        normalization_element = max(num_copies_a, num_copies_b)
            
        if shared_domain[:7] in anchor_domains:
            S_anchor += normalization_element
            domain_difference_anchor += sum_seq_dist
        else:
            S += normalization_element
            domain_difference += sum_seq_dist
            
    if S_anchor != 0 and S != 0:
        DSS_non_anchor = domain_difference / S
        DSS_anchor = domain_difference_anchor / S_anchor
        
        # Calculate proper, proportional weight to each kind of domain
        non_anchor_prct = S / (S + S_anchor)
        anchor_prct = S_anchor / (S + S_anchor)
        
        # boost anchor subcomponent and re-normalize
        non_anchor_weight = non_anchor_prct / (anchor_prct*anchorboost + non_anchor_prct)
        anchor_weight = anchor_prct*anchorboost / (anchor_prct*anchorboost + non_anchor_prct)

        # Use anchorboost parameter to boost percieved rDSS_anchor
        DSS = (non_anchor_weight*DSS_non_anchor) + (anchor_weight*DSS_anchor)
        
    elif S_anchor == 0:
        DSS_non_anchor = domain_difference / S
        DSS_anchor = 0.0
        
        DSS = DSS_non_anchor
        
    else: #only anchor domains were found
        DSS_non_anchor = 0.0
        DSS_anchor = domain_difference_anchor / S_anchor
        
        DSS = DSS_anchor
 
    DSS = 1-DSS #transform into similarity
 

    # ADJACENCY INDEX
    # calculates the Tanimoto similarity of pairs of adjacent domains
    
    if len(A_domlist[domA_start:domA_end]) < 2 or len(B_domlist[domB_start:domB_end]) < 2:
        AI = 0.0
    else:
        setA_pairs = set()
        for l in range(domA_start, domA_end-1):
            setA_pairs.add(tuple(sorted([A_domlist[l],A_domlist[l+1]])))
        
        setB_pairs = set()
        for l in range(domB_start, domB_end-1):
            setB_pairs.add(tuple(sorted([B_domlist[l],B_domlist[l+1]])))

        # same treatment as in Jaccard
        AI = len(setA_pairs & setB_pairs) / len(setA_pairs | setB_pairs)

    Distance = 1 - (Jaccardw * Jaccard) - (DSSw * DSS) - (AIw * AI)
    
    # This could happen due to numerical innacuracies
    if Distance < 0.0:
        if Distance < -0.000001: # this definitely is something else...
            print("Negative distance detected!")
            print(Distance)
            print("{} - {}".format(A, B))
            print("J: {}\tDSS: {}\tAI: {}".format(str(Jaccard), str(DSS), str(AI)))
            print("Jw: {}\tDSSw: {}\tAIw: {}".format(str(Jaccardw), str(DSSw), str(AIw)))
        Distance = 0.0
        
    rev = 0.0
    if reverse:
        rev = 1.0
        
    return Distance, Jaccard, DSS, AI, DSS_non_anchor, DSS_anchor, S, S_anchor, lcsStartA, lcsStartB, seedLength, rev


def launch_hmmalign(cores, domain_sequence_list):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    pool = Pool(cores, maxtasksperchild=32)
    pool.map(run_hmmalign, domain_sequence_list)
    pool.close()
    pool.join()
    
def run_hmmalign(domain_file):
    #domain_file already contains the full path, with the file extension
    domain_base = domain_file.split(os.sep)[-1][:-6]
    hmmfetch_pars = ["hmmfetch", os.path.join(pfam_dir,"Pfam-A.hmm.h3m"), domain_base]
    proc_hmmfetch = subprocess.Popen(hmmfetch_pars, stdout=subprocess.PIPE, shell=False)
    
    hmmalign_pars = ["hmmalign", "-o", domain_file.replace(".fasta",".stk"), "-", domain_file]
    proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout, stdout=subprocess.PIPE, shell=False)
    
    proc_hmmfetch.stdout.close()
    proc_hmmalign.communicate()[0]
    proc_hmmfetch.wait()
    
    if verbose:
        print(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))
    
    SeqIO.convert(domain_file[:-6]+".stk", "stockholm", domain_file[:-6]+".algn", "fasta")
    

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

def parseHmmScan(hmmscanResults, pfd_folder, pfs_folder, overlapCutoff):
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
            
    else:
        sys.exit("Error: hmmscan file " + outputbase + " was not found! (parseHmmScan)")

    return("")

def clusterJsonBatch(bgcs, pathBase, className, matrix, pos_alignments, 
                     cutoffs=[1.0], damping=0.8, clusterClans=False, 
                     clanCutoff=(0.5,0.8), htmlFolder=None):
    """BGC Family calling
    Uses csr sparse matrices to call Gene Cluster Families (GCFs) using Affinity
    Propagation.
    
    Cutoff values are in distance (i.e. if two clusters are further than cutoff 
    value, similarity is 0)
    Larger cutoff values are more permissive
    
    bgcs: ordered list of integers (ascending, unique, but not necessarily 
        consecutive) representing the index in the main clusterNames list. Every
        number here is an "external" index
    matrix: list of lists of idx0, idx1, d where the first two elements correspond
        to indices from `bgcs`
    pathBase: folder where GCF files will be deposited
    """
    numBGCs = len(bgcs)
    
    simDict = {} # dictionary of dictionaries
    # Doing this so it only has to go through the matrix once
    for row in matrix:
        gc1, gc2, distance = row
        
        if distance < 1.0:
            similarity = 1 - distance
        else:
            similarity = 0
        gcSimilarities = simDict.setdefault(gc1, {})
        gcSimilarities[gc2] = similarity

    clanClassificationCutoff, clanDistanceCutoff = clanCutoff
    if clusterClans and verbose:
        print('Clustering Clans Enabled with parameters clanClassificationCutoff: {}, clanDistanceCutoff: {}'.format(clanClassificationCutoff,clanDistanceCutoff))
    
    # External index number to Internal, consecutive, index
    bgcExt2Int = dict(zip(bgcs, range(numBGCs)))
    
    # Get info on all BGCs to export to .js for visualization
    bs_data = []
    bgcJsonDict = {}
    for bgc in bgcs:
        bgcName = clusterNames[bgc]
        bgcJsonDict[bgcName] = {}
        bgcJsonDict[bgcName]["id"] = bgcName
        bgcJsonDict[bgcName]["desc"] = bgc_info[bgcName].description
        bgcJsonDict[bgcName]["start"] = 1
        bgcJsonDict[bgcName]["end"] = bgc_info[bgcName].max_width
        
        pfdFile = os.path.join(pfd_folder, bgcName + ".pfd")
        fastaFile = os.path.join(bgc_fasta_folder, bgcName + ".fasta")
        
        orfDict = defaultdict(dict)
        
        ## read fasta file first to get orfs
        # We cannot get all the info exclusively from the pfd because that only
        # contains ORFs with predicted domains (and we need to draw empty genes
        # as well)
        for line in open(fastaFile):
            if line[0] == ">":
                header = line.strip()[1:].split(':')
                orf = header[0]
                if header[2]:
                    orfDict[orf]["id"] = header[2]
                elif header[4]:
                    orfDict[orf]["id"] = header[4]
                else:
                    orfDict[orf]["id"] = orf
                    
                ## broken gene goes into cluster, need this so js doesn't throw an error
                if int(header[6]) <= 1:
                    orfDict[orf]["start"] = 1
                else:
                    orfDict[orf]["start"] = int(header[6])
                    
                orfDict[orf]["end"] = int(header[7])
                
                if header[-1] == '+':
                    orfDict[orf]["strand"] = 1
                else:
                    orfDict[orf]["strand"] = -1
                    
                orfDict[orf]["domains"] = []
    
        ## now read pfd file to add the domains to each of the orfs
        for line in open(pfdFile):
            entry = line.split('\t')
            orf = entry[-1].strip().split(':')[0]
            pfamID = entry[5].split('.')[0]
            orfDict[orf]["domains"].append({'code': pfamID, 'start': int(entry[3]), 'end': int(entry[4]), 'bitscore': float(entry[1])})
        # order of ORFs is important here because I use it to get a translation
        # between the "list of ORFs with domains" to "list of all ORFs" later on
        bgcJsonDict[bgcName]['orfs'] = sorted(orfDict.values(), key=itemgetter("start"))
    bs_data = [bgcJsonDict[clusterNames[bgc]] for bgc in bgcs]
    
    ## Write html output folder structure (and update bigscape_results.js) for this module
    assert os.path.isdir(htmlFolder)
    module_html_path = os.path.join(htmlFolder, className)
    create_directory(module_html_path, "Network HTML", False)
    with open(os.path.join(module_html_path, "bs_data.js"), "w") as bs_data_js:
        bs_data_js.write("var bs_data={};\n".format(json.dumps(bs_data, indent=4, separators=(',', ':'), sort_keys=True)))
    shutil.copy(os.path.join(os.path.realpath(os.path.dirname(__file__)), "html_template", "index_html"), os.path.join(module_html_path, "index.html"))
    
    for cutoff in cutoffs:
        #simMatrix = lil_matrix((numBGCs, numBGCs), dtype=np.float32)
        simMatrix = np.zeros((numBGCs, numBGCs), dtype=np.float32)
        for bgc1 in bgcs:
            # first make sure it is similar to itself
            simMatrix[bgcExt2Int[bgc1],bgcExt2Int[bgc1]] = 1
            for bgc2 in simDict.get(bgc1,{}).keys():
                # you might get 0 values if there were matrix entries under the
                # cutoff. No need to input these in the sparse matrix
                if simDict[bgc1][bgc2] > 1-cutoff:
                    # Ensure symmetry
                    simMatrix[bgcExt2Int[bgc1], bgcExt2Int[bgc2]] = simDict[bgc1][bgc2]
                    simMatrix[bgcExt2Int[bgc2], bgcExt2Int[bgc1]] = simDict[bgc1][bgc2]
        if verbose:
            print("  Clustering (c=" + str(cutoff) + ")")
        #labels = pysapc.SAP(damping=damping, max_iter=500,
                            #preference='min').fit_predict(simMatrix)
        af = AffinityPropagation(damping=damping, max_iter=500).fit(simMatrix)
        # labels_: labels for each node. Same label means same cluster. Value is also
        #   the index in the cluster_centers_indices_ list. Goes from 0 to number of 
        #   clusters-1
        # cluster_centers_indices_: list of exemplars
        labels = [af.cluster_centers_indices_[af.labels_[i]] for i in range(numBGCs)]
        
        if verbose:
            print("   ...done")

        bs_distances = [[float("{:.3f}".format(simMatrix[row, col])) for col in 
                         range(row+1)] for row in range(numBGCs)]
        
        
        
        familiesDict = defaultdict(list)
        for i in range(numBGCs):
            familiesDict[bgcs[labels[i]]].append(bgcs[i])
        familyIdx = sorted(familiesDict.keys()) # identifiers for each family
        
        #write_error = False
        #for ex in familiesDict:
            #if ex not in familiesDict[ex]:
                #write_error = True
                #print(ex, familiesDict[ex])
        #if write_error:
            #print("Errors where found in scikit-learn, writing errorfile")
            #with open("errorfile.txt","w") as err:
                #err.write("Original bgc index\tInternal indexes (0-numBGCs-1)\tLabel (internal index)\n")
                #for i in range(numBGCs):
                    #err.write("{}\t{}\t{}\n".format(bgcs[i],i,labels[i]))
            #sys.exit("stop...")

        ##
        ## Get conserved domain core information to build phylogenetic tree
        ## 
        gcf_trees_path = os.path.join(pathBase, "GCF_trees")
        if not os.path.exists(gcf_trees_path):
            os.makedirs(gcf_trees_path)
            
        newick_trees = {} # key:original bgc index
        for exemplar_idx in familiesDict:
            exemplar = clusterNames[exemplar_idx]
            gcf = familiesDict[exemplar_idx][:]
            if len(gcf) < 3:
                newick_trees[exemplar_idx] = "({});".format(",".join([str(bgcExt2Int[x])+":0.0" for x in gcf]))
                
                #print(newick_trees[exemplar_idx])
                #TODO make some default alignment data to send to the json file
                continue
            
            domain_sets = {}
            
            # make a frequency table (not counting copies):
            frequency_table = defaultdict(int)
            for bgc in gcf:
                domain_sets[bgc] = set(DomainList[clusterNames[bgc]])
                for domain in domain_sets[bgc]:
                    frequency_table[domain] += 1
            
            # Remove all PKS/NRPS domains
            # TODO but this was intended to be done when we were considering
            # directly using Corason. Should these domains still be removed?
            # We are now able to include them accurately with the Hungarian
            # matching
            # If implemented, fill nrps_pks_domains first!
            #domains_in_gcf = set(frequency_table.keys())
            #for erase_domain in (nrps_pks_domains & domains_in_gcf):
                #del frequency_table[erase_domain]
                
            # Find the set of [(tree domains)]. They should 1) be in the exemplar
            # and 2) appear with the most frequency. Iterate over the different
            # frequencies (descending) until set is not empty
            tree_domains = set()
            frequencies = sorted(set(frequency_table.values()), reverse=True)            
            
            # first try with domain(s) with max frequency, even if it's just one
            f = 0 
            while len(tree_domains) == 0 and f < len(frequencies):
                for domain in frequency_table:
                    if frequency_table[domain] == frequencies[f] and domain in domain_sets[exemplar_idx]:
                            tree_domains.add(domain)
                    
                f += 1
            
            
            if len(tree_domains) == 1:
                print("  Warning: core shared domains for GCF {} consists of a single domain ({})".format(exemplar_idx, [x for x in tree_domains][0]))
            
            # Get the alignments of the core domains
            alignments = {}
            # initialize every sequence alignment entry. Don't do defaultdict!
            alignments[exemplar_idx] = ""
            
            out_of_tree_bgcs = [] # bgcs that don't share a common domain core
            delete_list = []
            gcf.remove(exemplar_idx) # separate exemplar from the rest of the bgcs
            for bgc in gcf:
                alignments[bgc] = ""

            match_dict = {}
            for domain in tree_domains:
                #print(domain)
                specific_domain_list_A = BGCs[exemplar][domain]
                num_copies_a = len(specific_domain_list_A)
                for exemplar_domain_copy in specific_domain_list_A:
                    alignments[exemplar_idx] += AlignedDomainSequences[exemplar_domain_copy]
                
                for bgc in gcf:
                    match_dict.clear()
                    if domain not in domain_sets[bgc]:
                        print("BGC {} does not share a common domain core (domain: {})".format(clusterNames[bgc], domain))
                        out_of_tree_bgcs.append(bgc)
                        delete_list.append(bgc)
                    else:
                        specific_domain_list_B = BGCs[clusterNames[bgc]][domain]
                        
                        num_copies_b = len(specific_domain_list_B)
                        
                        DistanceMatrix = np.ndarray((num_copies_a,num_copies_b))
                        for domsa in range(num_copies_a):
                            for domsb in range(num_copies_b):
                                # TODO NOT taking into consideration any LCS slicing
                                # i.e. we're comparing ALL copies of this domain
                                sequence_tag_a = specific_domain_list_A[domsa]
                                sequence_tag_b = specific_domain_list_B[domsb]
                                
                                aligned_seqA = AlignedDomainSequences[sequence_tag_a]
                                aligned_seqB = AlignedDomainSequences[sequence_tag_b]
                                
                                seq_length = 0
                                matches = 0
                                gaps = 0
                                
                                seq_length = len(aligned_seqA)
                                
                                for position in range(seq_length):
                                    if aligned_seqA[position] == aligned_seqB[position]:
                                        if aligned_seqA[position] != "-":
                                            matches += 1
                                        else:
                                            gaps += 1
                                            
                                DistanceMatrix[domsa][domsb] = 1 - ( matches/(seq_length-gaps) )

                        BestIndexes = linear_sum_assignment(DistanceMatrix)
                        # at this point is not ensured that we have the same order
                        # for the exemplar's copies (rows in BestIndexes)
                        # ideally they should go from 0-numcopies. Better make sure
                        for x in range(len(BestIndexes[0])):
                            match_dict[BestIndexes[0][x]] = BestIndexes[1][x]
                            
                        for copy in range(num_copies_a):
                            try:
                                alignments[bgc] += AlignedDomainSequences[specific_domain_list_B[match_dict[copy]]]
                            except KeyError:
                                # This means that this copy of exemplar did not
                                # have a match in bgc (i.e. bgc has less copies
                                # of this domain than exemplar)
                                alignments[bgc] += "-"*seq_length
                        
                for bgc in delete_list:
                    try:
                        del alignments[bgc]
                    except KeyError:
                        # if using defaultdict, the except will always add the
                        # key when except'ing
                        pass
                del delete_list[:]
                 
            # need this to change the labels in the trees that are read from files
            bgc_name_to_idx = {}
            
            # save compiled alignments of the GCF domain core as fastas
            alignment_file_path = os.path.join(gcf_trees_path,"GCF_c{:4.2f}_{:05d}_alignment.fasta".format(cutoff,exemplar_idx))
            with open(alignment_file_path, "w") as gcf_alignment_file:
                gcf_alignment_file.write(">{}\n{}\n".format(exemplar, alignments[exemplar_idx]))
                bgc_name_to_idx[exemplar] = exemplar_idx
                
                for bgc in alignments:
                    if bgc != exemplar_idx:
                        gcf_alignment_file.write(">{}\n{}\n".format(clusterNames[bgc], alignments[bgc]))
                        bgc_name_to_idx[clusterNames[bgc]] = bgc
            

            # launch fasttree to make tree
            if verbose:
                print("  Working GCF {}, cutoff {}".format(exemplar_idx, cutoff))
            newick_file_path = os.path.join(gcf_trees_path, "GCF_c{:4.2f}_{:05d}.newick".format(cutoff,exemplar_idx))
            with open(newick_file_path, "w") as newick_file:
                command = ["fasttree", "-nopr", "-quiet", alignment_file_path]
                p = subprocess.Popen(command, stdout=newick_file, shell=False)
                p.wait() # only with process has terminated will the file be ready

            
            if not os.path.isfile(newick_file_path):
                print(newick_file_path)
                sys.exit(" ERROR: newick file not created (GCF_c{:4.2f}_{:05d})".format(cutoff,exemplar_idx))
            else:
                # change bgc names to (internal) indices
                #with open(newick_file_path,"r") as newick_file:
                    #newick = newick_file.read().strip()
                    #print(newick)
                    #print("")
                    #for name in bgc_name_to_idx:
                        #newick = newick.replace(name, str(bgcExt2Int[bgc_name_to_idx[name]]))
                    #newick_trees[exemplar_idx] = newick
                    
                    #print(newick)
                    #print("")
                    
                with open(newick_file_path,"r") as newick_file:
                    tree = Phylo.read(newick_file, 'newick')
                    #print(tree.format("newick"))
                    #print("")
                    tree.root_at_midpoint()
                    newick = tree.format("newick")
                    for name in bgc_name_to_idx:
                        #print("Replace {} with {}".format(name, str(bgcExt2Int[bgc_name_to_idx[name]])))
                        newick = newick.replace(name, str(bgcExt2Int[bgc_name_to_idx[name]]))
                    #print(newick)
                    newick_trees[exemplar_idx] = newick
       
        ### Use the 0.5 distance cutoff to cluster clans by default
        if clusterClans and cutoff == clanClassificationCutoff:
            # Detect if there's only 1 GCF. It makes pySAPC crash
            if len(familyIdx) == 1:
                #print("  Gene Cluster Clans: Only 1 GCF")
                clanLabels = [1]
                continue
            
            famSimMatrix = lil_matrix((len(familyIdx), len(familyIdx)), dtype=np.float32)
            #famSimMatrix = np.zeros((len(familyIdx), len(familyIdx)), dtype=np.float32)
            familiesExt2Int = {gcfExtIdx:gcfIntIdx for gcfIntIdx,gcfExtIdx in enumerate(familyIdx)}
                                                                                        
            for familyI, familyJ in [tuple(sorted(combo)) for combo in combinations(familyIdx, 2)]:
                famSimilarities = []
                # currently uses the average distance of all average distances
                # between bgc from gcf I to all bgcs from gcf J
                for bgcI in familiesDict[familyI]:
                    similarities = [simMatrix[bgcExt2Int[bgcI], bgcExt2Int[bgcJ]] for bgcJ in familiesDict[familyJ]]
                    famSimilarities.append(sum(similarities, 0.0) / len(similarities))
                try:
                    familySimilarityIJ = sum(famSimilarities, 0.0)/len(famSimilarities)
                except ZeroDivisionError:
                    familySimilarityIJ = 0.0
                if familySimilarityIJ > 1 - clanDistanceCutoff:
                    # Ensure symmetry
                    famSimMatrix[familiesExt2Int[familyI], familiesExt2Int[familyJ]] = familySimilarityIJ
                    famSimMatrix[familiesExt2Int[familyJ], familiesExt2Int[familyI]] = familySimilarityIJ
            # add main diagonal
            for family in range(familyIdx):
                famSimMatrix[family,family] = 1.0
            
            clanLabels = pysapc.SAP(damping=damping, max_iter=500,
                                preference='min').fit_predict(famSimMatrix)
        else:
            clanLabels = []
            
        if len(clanLabels) > 0:
            clansDict = defaultdict(list)
            for i in range(familyIdx):
                clansDict[clanLabels[i]].append(familyIdx[i])

            fam2clan = dict(zip(familyIdx,clanLabels))
        
            bs_families = [{"id": "FAM_{:05d}".format(family), 
                            'members': [bgcExt2Int[member] for member in members], 
                            "clan": "CLAN_{:03d}".format(fam2clan[family])}
                           for family, members in familiesDict.items()]
            bs_clans = [{"id": "CLAN_{:03d}".format(clan), 'members': members}
                           for clan, members in clansDict.items()]
        else:
            bs_families = [{"id": "FAM_{:05d}".format(family), 
                            'members': [bgcExt2Int[member] for member in members], }
                           for family, members in familiesDict.items()]
        
        # Positional alignment information is based on DomainCountGene, which
        # does not contain empty genes (i.e. with no domains). 
        domainGenes2allGenes = {}
        
        ## BGC Family alignment information
        #print(pos_alignments)
        bs_families_alignment = []
        for family, members in familiesDict.items():
            #print("FAMILY {}".format(family))
            for bgc in members:
                #print(bgc, clusterNames[bgc])
                domainGenes2allGenes[bgc] = {}
                has_domains = 0
                for orf in range(len(bs_data[bgcExt2Int[bgc]]["orfs"])):
                    #print("{} domains in orf {}".format(len(bs_data[bgcExt2Int[bgc]]["orfs"][orf]["domains"]), orf))
                    if len(bs_data[bgcExt2Int[bgc]]["orfs"][orf]["domains"]) > 0:
                        domainGenes2allGenes[bgc][has_domains] = orf
                        has_domains += 1
                #print(len(bs_data[bgcExt2Int[bgc]]["orfs"]))
            #print(domainGenes2allGenes)
            assert (len(members) > 0), "Error: bs_families[{}] have no members, something went wrong?".format(fam_idx)
            
            ref_genes_ = set()
            aln = []
            for bgc in members:
                if bgc == family:
                    aln.append([ [gene_num, 0] for gene_num in range(len(bs_data[bgcExt2Int[family]]["orfs"]))])
                else:
                    try:
                        a, b, length, reverse = pos_alignments[family][bgc]
                    except:
                        b, a, length, reverse = pos_alignments[bgc][family]
                        
                        if reverse:
                            # special case. bgc was reference (first) in lcs
                            a = domainGenes2allGenes[family][len(DomainCountGene[clusterNames[family]])-a-length]
                            b = domainGenes2allGenes[bgc][b+length-1] # -1 go to 0-index
                            
                            #if clusterNames[family] == "P207.1.cluster032" and clusterNames[bgc] == "LGCW01000001.1.cluster034":
                                #print(length)
                                #print(DomainCountGene[clusterNames[family]], len(DomainCountGene[clusterNames[family]]))
                                #print(domainGenes2allGenes[family])
                                #print(DomainCountGene[clusterNames[bgc]], len(DomainCountGene[clusterNames[bgc]]))
                                #print(domainGenes2allGenes[bgc])
                        else:
                            a = domainGenes2allGenes[family][a]
                            b = domainGenes2allGenes[bgc][b]
                    else:
                        #print(clusterNames[family], a, clusterNames[bgc], b)
                        a = domainGenes2allGenes[family][a]
                        
                        if reverse:
                            b = domainGenes2allGenes[bgc][len(DomainCountGene[clusterNames[bgc]])-b-1]
                            #print(a, b)
                        else:
                            b = domainGenes2allGenes[bgc][b]
                            
                            
                    ref_genes_.add(a)
                    bgc_algn = []
                    
                    for gene_num in range(len(bs_data[bgcExt2Int[bgc]]["orfs"])):
                        if gene_num == b: # this is the reference gene for this bgc
                            if reverse:
                                bgc_algn.append([a, -100])
                            else:
                                bgc_algn.append([a, 100])
                        else:
                            bgc_algn.append([-1, 100])
                    
                    aln.append(bgc_algn)
                        
            ref_genes = list(ref_genes_)
            
            fam_alignment = {
                "id" : "FAM_{:05d}".format(family),
                "ref" : bgcExt2Int[family],
                "newick" : newick_trees[family],
                "ref_genes" : ref_genes,
                "aln" : aln
            }
            bs_families_alignment.append(fam_alignment)
        ## End of BGC Family alignment information

        # column1: BGC, column2: clustering pseudo family
        if verbose:
            print("  Writing clustering file")
        clustering_file_path = os.path.join(pathBase, "{}_clustering_c{:4.2f}.tsv".format(className, cutoff))
        with open(clustering_file_path, "w") as clustering_file:
            clustering_file.write('#BGC Name\tFamily Number\n')
            for family in familyIdx:
                for bgc in familiesDict[family]:
                    clustering_file.write("{}\t{}\n".format(clusterNames[bgc], family))

        ## Write bgc_networks.js		
        with open(os.path.join(module_html_path, "bs_networks.js"), "w") as bs_networks_js:
            bs_networks_js.write("var bs_similarity={};\n".format(json.dumps(bs_distances, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families={};\n".format(json.dumps(bs_families, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families_alignment={};\n".format(json.dumps(bs_families_alignment, indent=4, separators=(',', ':'), sort_keys=True)))
            if len(clanLabels) > 0:
                bs_networks_js.write("var bs_clans={};\n".format(json.dumps(bs_clans, indent=4, separators=(',', ':'), sort_keys=True)))
        
        if len(clanLabels) > 0:
            if verbose:
                print("   Writing Clans file")
            clans_file_path = os.path.join(pathBase, "{}_clans_{:4.2f}_{:4.2f}.tsv".format(className,clanClassificationCutoff,clanDistanceCutoff))
            with open(clans_file_path,'w') as clansFile:
                clansFile.write('#BGC Name\tClan Number\tFamily Number\n')
                for clan in clansDict.keys():
                    for family in clansDict[clan]:
                        for bgc in familiesDict[family]:
                            clansFile.write("{}\t{}\t{}\n".format(clusterNames[bgc], clan, family))
                    
    return


class FloatRange(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end
    def __repr__(self):
        return '{}-{}'.format(self.start, self.end)

def CMD_parser():
    parser = ArgumentParser(prog="BiG-SCAPE")
    
    parser.add_argument("-i", "--inputdir", dest="inputdir", 
                        default=os.path.dirname(os.path.realpath(__file__)),
                        help="Input directory of gbk files, if left empty, all \
                        gbk files in current and lower directories will be used.")

    parser.add_argument("-o", "--outputdir", dest="outputdir", default="", 
                        required=True, help="Output directory, this will contain \
                        all output data files.")
    
    parser.add_argument("--pfam_dir", dest="pfam_dir",
                      default=os.path.dirname(os.path.realpath(__file__)), 
                      help="Location of hmmpress-processed Pfam files. Default \
                      is same location of BiG-SCAPE")
    
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(),
                      help="Set the number of cores the script may use (default:\
                      use all available cores)")
    
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", 
                        default=False, help="Prints more detailed information. \
                        Toggle to activate.")
    
    parser.add_argument("--include_singletons", dest="include_singletons", 
                        action="store_true", default=False, help="Include nodes \
                        that have no edges to other nodes from the network. \
                        Toggle to activate.")
    
    parser.add_argument("-d", "--domain_overlap_cutoff", 
                        dest="domain_overlap_cutoff", default=0.1, help="Specify\
                        at which overlap percentage domains are considered to \
                        overlap. Domain with the best score is kept \
                        (default=0.1).")
    
    parser.add_argument("-m", "--min_bgc_size", dest="min_bgc_size", default=0,
                      help="Provide the minimum size of a BGC to be included in\
                      the analysis. Default is 0 base pairs")
    
    parser.add_argument("-s", "--samples", dest="samples", action="store_true", 
                        default=False, help="Separate the input files into \
                        samples according to their containing folder within the \
                        input folder. Toggle to activate")
    
    parser.add_argument("--no_all", dest="no_all", action="store_true", 
                        default=False, help="By default, BiG-SCAPE uses a single\
                        data set comprised of all input files available \
                        recursively within the input folder. Toggle to \
                        disactivate this behaviour (in that case, if the \
                        --samples parameter is not activated, BiG-SCAPE will not\
                        create any network file)")
    
    parser.add_argument("--mix", dest="mix", action="store_true", default=False, 
                        help="By default, BiG-SCAPE separates the analysis \
                        according to the BGC product (PKS Type I, NRPS, RiPPs, etc.) \
                        and will create network directories for each class. \
                        Toggle to include an analysis mixing all classes")
    
    parser.add_argument("--no_classify", dest="no_classify", action="store_true", 
                        default=False, help="By default, BiG-SCAPE classifies \
                        the output files analysis based on the BGC product. \
                        Toggle to deactivate (note that if the --mix parameter \
                        is not activated, BiG-SCAPE will not create any network \
                        file).")
    
    parser.add_argument("--banned_classes", nargs='+', dest="banned_classes", 
                        default=[], choices=["PKSI", "PKSother", "NRPS", "RiPPs", 
                                             "Saccharides", "Terpene", 
                                             "PKS-NRP_Hybrids", "Others"], 
                        help="Classes that should NOT be included in the \
                        classification. E.g. \"--banned_classes PKSI PKSOther\"")

    
    parser.add_argument("--clans", dest="clans",action="store_true", 
                        default=False, help="BiG-SCAPE will perform a second \
                        layer of clustering and attempt to group families \
                        assigned from clustering with cutoff of 0.5 to clans")

    parser.add_argument("--clan_cutoff",dest="clan_cutoff",default=[0.5,0.8], 
                        type=float, choices=[FloatRange(0.0, 1.0)],nargs=2,
                        help="Cutoff Parameters for which clustering families \
                        into clans will be performed in raw distance. First \
                        value is the cutoff value family assignments for BGCs \
                        used in clan clustering (default: 0.5). Second value is \
                        the cutoff value for clustering families into clans \
                        (default: 0.8). Average linkage for BGCs in a family is\
                        used for distances between families. Example: \
                        --clan_cutoff 0.5 0.8)")

    parser.add_argument("--hybrids", dest="hybrids", action="store_true", 
                        default=False, help="Toggle to also add BGCs with hybrid\
                        predicted products from the PKS/NRPS Hybrids and Others\
                        classes to each subclass (e.g. a 'terpene-nrps' BGC from\
                        Others would be added to the Terpene and NRPS classes")
    
    parser.add_argument("--mode", dest="mode", default="global", choices=["global",
                            "lcs", "auto"], help="Alignment mode for each pair of\
                            gene clusters. 'global' (default) the whole list of\
                            domains of each BGC are compared; 'lcs': Longest\
                            Common Subcluster mode. Redefine the subset of the \
                            domains used to calculate distance by trying to find\
                            the longest slice of common domain content per gene\
                            in both BGCs, then expand each slice.\
                            'auto' use LCS when at least one of the BGCs in each\
                            pair has the 'contig_edge' annotation from antiSMASH\
                            v4+, otherwise use global mode on that pair")
    
    parser.add_argument("--anchorfile", dest="anchorfile", 
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),"anchor_domains.txt"),
                        help="Provide a custom location for the anchor domains \
                        file, default is anchor_domains.txt.")
    
    parser.add_argument("--exclude_gbk_str", dest="exclude_gbk_str", 
                        default="final", nargs="+",
                        help="If any string in this list occurs in the gbk \
                        filename, this file will not be used for the analysis.\
                        (default: final)")
                    
    parser.add_argument("--force_hmmscan", dest="force_hmmscan", action="store_true", 
                        default=False, help="Force domain prediction using \
                        hmmscan even if BiG-SCAPE finds processed domtable files\
                        (e.g. to use a new version of PFAM).")
    parser.add_argument("--skip_ma", dest="skip_ma", action="store_true", 
                        default=False, help="Skip multiple alignment of domains'\
                        sequences. Use if alignments have been generated in a \
                        previous run.")
    parser.add_argument("--skip_all", dest="skip_all", action="store_true",
                      default = False, help ="Only generate new network files.")
    parser.add_argument("--cutoffs", dest="cutoffs", nargs="+", default=[0.10, 
                    0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 
                    0.65, 0.70, 0.75, 0.80], type=float, choices=[FloatRange(0.0, 
                    1.0)], help="Generate networks using multiple raw distance \
                    cutoff values, example: --cutoffs 0.1, 0.25, 0.5, 1.0. Default: \
                    all values from 0.10 to 0.80 with 0.05 intervals.")

    parser.add_argument("--version", action="version", version="%(prog)s 201709")

    return parser.parse_args()


if __name__=="__main__":
    options = CMD_parser()
    
    class bgc_data:
        def __init__(self, accession_id, description, product, records, max_width, organism, taxonomy, biosynthetic_genes, contig_edge):
            # These two properties come from the genbank file:
            self.accession_id = accession_id
            self.description = description
            # AntiSMASH predicted class of compound:
            self.product = product
            # number of records in the genbank file (think of multi-locus BGCs):
            self.records = records
            # length of largest record (it will be used for ArrowerSVG):
            self.max_width = int(max_width)
            # organism
            self.organism = organism
            # taxonomy as a string (of comma-separated values)
            self.taxonomy = taxonomy
            # Internal set of tags corresponding to genes that AntiSMASH marked 
            # as "Kind: Biosynthetic". It is formed as
            # clusterName + "_ORF" + cds_number + ":gid:" + gene_id + ":pid:" + protein_id + ":loc:" + gene_start + ":" + gene_end + ":strand:" + {+,-}
            self.biosynthetic_genes = biosynthetic_genes
            # AntiSMASH 4+ marks BGCs that sit on the edge of a contig
            self.contig_edge = contig_edge

    
    if options.outputdir == "":
        print("please provide a name for an output folder using parameter -o or --outputdir")
        sys.exit(0)
    
    global anchor_domains
    if os.path.isfile(options.anchorfile):
        anchor_domains = get_anchor_domains(options.anchorfile)
    else:
        print("File with list of anchor domains not found")
        anchor_domains = set()
    
    global bgc_class_weight
    global AlignedDomainSequences
    global DomainList
    global DomainCountGene
    global corebiosynthetic_position
    global verbose
    global BGCs
    
    # contains the type of the final product of the BGC (as predicted by AntiSMASH), 
    # as well as the definition line from the BGC file. Used in the final network files.
    global output_folder

    global pfam_dir
    global timings_file
    global cores
    global mode
    
    global run_name

    global clusterNames, bgcClassNames
    
    include_singletons = options.include_singletons
    
    cores = int(options.cores)
    
    options_all = not options.no_all
    options_samples = options.samples
    
    options_mix = options.mix
    options_classify = not options.no_classify
    
    mode = options.mode
    
    cutoff_list = options.cutoffs
    for c in cutoff_list:
        if c <= 0.0 or c > 1.0:
            print(" Removing invalid cutoff value {}".format(str(c)))
            cutoff_list.remove(c)
            
            
    # if we want to classify by clans make sure that the clanCutoff is included in the cutoffs to do AP in
    if options.clans:
        fc, cc = options.clan_cutoff
        if c not in cutoff_list:
            if fc <= 0.0 or fc > 1.0:
                sys.exit("Error: invalid cutoff value for GCF calling")
            else:
                cutoff_list = sorted(cutoff_list.append(fc))
            
        if cc <= 0.0 or fc > 1.0:
            sys.exit("Error: invalid cutoff value for GCC calling")
        

    output_folder = str(options.outputdir)
    
    pfam_dir = str(options.pfam_dir)
    h3f = os.path.join(pfam_dir, "Pfam-A.hmm.h3f")
    h3i = os.path.join(pfam_dir, "Pfam-A.hmm.h3i")
    h3m = os.path.join(pfam_dir, "Pfam-A.hmm.h3m")
    h3p = os.path.join(pfam_dir, "Pfam-A.hmm.h3p")
    if not (os.path.isfile(h3f) and os.path.isfile(h3i) and os.path.isfile(h3m) and os.path.isfile(h3p)):
        print("One or more of the necessary Pfam files (.h3f, .h3i, .h3m, .h3p) were not found")
        if os.path.isfile(os.path.join(pfam_dir, "Pfam-A.hmm")):
            print("Please use hmmpress with Pfam-A.hmm")
        else:
            print("Please download the latest Pfam-A.hmm file from http://pfam.xfam.org/")
            print("Then use hmmpress on it, and use the --pfam_dir parameter to point to the location of the files")
        sys.exit()

                    
    verbose = options.verbose
    
    run_mode_string = ""
    networks_folder_all = "networks_all"
    networks_folder_samples = "networks_samples"    
    if options.hybrids:
        networks_folder_all += "_hybrids"
        networks_folder_samples += "_hybrids"
        run_mode_string += "_hybrids"
    if mode == "auto":
        networks_folder_all += "_auto"
        networks_folder_samples += "_auto"
        run_mode_string += "_auto"
    elif mode == "lcs":
        networks_folder_all += "_lcs"
        networks_folder_samples += "_lcs"
        run_mode_string += "_glocal"
    else:
        run_mode_string += "_full"
    
    if options.skip_all and options.skip_ma:
        print("Overriding --skip_ma with --skip_all parameter")
        options.skip_hmmscan = False
        options.skip_ma = False

    time1 = time.time()

    run_name = "{}{}".format(time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime()), run_mode_string)
    
    # Make the following available for possibly deleting entries within parseHmmScan
    global genbankDict, gbk_files, sampleDict, clusters, baseNames
    
    
    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    print("\n\n   - - Processing input files - -")
    
    create_directory(output_folder, "Output", False)

    # logs
    log_folder = os.path.join(output_folder, "logs")
    create_directory(log_folder, "Logs", False)
    write_parameters(log_folder, sys.argv)

    # cached stuff
    cache_folder = os.path.join(output_folder, "cache")
    bgc_fasta_folder = os.path.join(cache_folder, "fasta")
    domtable_folder = os.path.join(cache_folder, "domtable")    
    pfs_folder = os.path.join(cache_folder, "pfs")
    pfd_folder = os.path.join(cache_folder, "pfd")    
    domains_folder = os.path.join(cache_folder, "domains")
    create_directory(cache_folder, "Cache", False)
    create_directory(bgc_fasta_folder, "BGC fastas", False)
    create_directory(domtable_folder, "Domtable", False)
    create_directory(domains_folder, "Domains", False)    
    create_directory(pfs_folder, "pfs", False)
    create_directory(pfd_folder, "pfd", False)

    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    bgc_info = {} # Stores, per BGC: predicted type, gbk Description, number of records, width of longest record, GenBank's accession, Biosynthetic Genes' ids
    genbankDict = get_gbk_files(options.inputdir, output_folder, bgc_fasta_folder, int(options.min_bgc_size), options.exclude_gbk_str, bgc_info)

    # clusters and sampleDict contain the necessary structure for all-vs-all and sample analysis
    clusters = genbankDict.keys()
    
    sampleDict = {} # {sampleName:set(bgc1,bgc2,...)}
    gbk_files = [] # raw list of gbk file locations
    for (cluster, (path, clusterSample)) in genbankDict.items():
        gbk_files.append(path)
        for sample in clusterSample:
            clustersInSample = sampleDict.get(sample, set())
            clustersInSample.add(cluster)
            sampleDict[sample] = clustersInSample

    
    print("\nCreating output directories")
    svg_folder = os.path.join(output_folder, "SVG")    
    create_directory(svg_folder, "SVG", False)
    network_folder = os.path.join(output_folder, "network_files")
    create_directory(network_folder, "Networks", False)

    print("\nTrying threading on {} cores".format(str(cores)))
    
    """BGCs -- 
    dictionary of this structure:
    BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""     
    BGCs = {} #will contain the BGCs
    
    
    # Weights in the format J, DSS, AI, anchorboost
    # Generated with optimization results 2016-12-05. 
    # Used the basic list of 4 anchor domains.
    bgc_class_weight = {}
    bgc_class_weight["PKSI"] = (0.22, 0.76, 0.02, 1.0)
    bgc_class_weight["PKSother"] = (0.0, 0.32, 0.68, 4.0)
    bgc_class_weight["NRPS"] = (0.0, 1.0, 0.0, 4.0)
    bgc_class_weight["RiPPs"] = (0.28, 0.71, 0.01, 1.0)
    bgc_class_weight["Saccharides"] = (0.0, 0.0, 1.0, 1.0)
    bgc_class_weight["Terpene"] = (0.2, 0.75, 0.05, 2.0)
    bgc_class_weight["PKS-NRP_Hybrids"] = (0.0, 0.78, 0.22, 1.0)
    bgc_class_weight["Others"] = (0.01, 0.97, 0.02, 4.0)
    
    #define which classes will be analyzed (if in the options_classify mode)
    valid_classes = set()
    for key in bgc_class_weight:
        valid_classes.add(key.lower())
    user_banned_classes = set([a.strip().lower() for a in options.banned_classes])
    valid_classes = valid_classes - user_banned_classes

    bgc_class_weight["mix"] = (0.2, 0.75, 0.05, 2.0) # default when not separating in classes
    BGC_classes = defaultdict(list)
    # mix class will always be the last element of the tuple
    bgcClassNames = tuple(sorted(list(bgc_class_weight)) + ["mix"])
    assert bgcClassNames[-1] == 'mix'

    bgcClassName2idx = dict(zip(bgcClassNames,range(len(bgcClassNames))))

    AlignedDomainSequences = {} # Key: specific domain sequence label. Item: aligned sequence
    DomainList = {} # Key: BGC. Item: ordered list of domains
    
    # Key: BGC. Item: ordered list of simple integers with the number of domains
    # in each gene
    # Instead of `DomainCountGene = defaultdict(list)`, let's try arrays of 
    # unsigned ints
    DomainCountGene = {}
    # list of gene-numbers that have a hit in the anchor domain list. Zero based
    corebiosynthetic_position = {}
    # list of +/- orientation 
    BGCGeneOrientation = {}
    
    # to avoid multiple alignment if there's only 1 seq. representing a particular domain
    sequences_per_domain = {}
    
    
    ### Step 2: Run hmmscan
    print("\nPredicting domains using hmmscan")
    
    baseNames = set(clusters)
    
    # All available fasta files (could be more than it should if reusing output folder)
    allFastaFiles = set(glob(os.path.join(bgc_fasta_folder,"*.fasta")))
    
    # fastaFiles: all the fasta files that should be there 
    # (i.e. correspond to the input files)
    fastaFiles = set()
    for name in baseNames:
        fastaFiles.add(os.path.join(bgc_fasta_folder, name+".fasta"))

    # fastaBases: the actual fasta files we have that correspond to the input
    fastaBases = allFastaFiles.intersection(fastaFiles)
    
    # Verify that all input files had their fasta sequences extracted
    if len(fastaFiles - fastaBases) > 0:
        sys.exit("Error! The following files did NOT have their fasta sequences extracted: " + ", ".join(fastaFiles - fastaBases))
    
    # Make a list of all fasta files that need to be processed
    # (i.e., they don't yet have a corresponding .domtable)
    if options.force_hmmscan:
        # process all files, regardless of whether they already existed
        task_set = fastaFiles
        print(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
    else:
        # find already processed files
        alreadyDone = set()
        for fasta in fastaFiles:
            outputbase  = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(domtable_folder,outputbase + '.domtable')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(fasta)
        task_set = fastaFiles - alreadyDone
        if len(task_set) == 0:
            print(" All fasta files had already been processed")
        elif len(alreadyDone) > 0:
            if len(task_set) < 20:
                print(" Warning! The following NEW fasta file(s) will be processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in task_set)))
            else:
                print(" Warning: {} NEW fasta files will be processed".format(str(len(task_set))))
        else:
            print(" Predicting domains for {} fasta files".format(str(len(fastaFiles))))
        
    pool = Pool(cores,maxtasksperchild=1)
    for fastaFile in task_set:
        pool.apply_async(runHmmScan,args=(fastaFile, pfam_dir, domtable_folder, verbose))
    pool.close()
    pool.join()
    print(" Finished generating domtable files.")

    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    print("\nParsing hmmscan domtable files")
    
    # All available domtable files
    allDomtableFiles = set(glob(os.path.join(domtable_folder,"*.domtable")))
    
    # domtableFiles: all domtable files corresponding to the input files
    domtableFiles = set()
    for name in baseNames:
        domtableFiles.add(os.path.join(domtable_folder, name+".domtable"))
    
    # domtableBases: the actual set of input files with coresponding domtable files
    domtableBases = allDomtableFiles.intersection(domtableFiles)
    
    # Verify that all input files have a corresponding domtable file
    if len(domtableFiles - domtableBases) > 0:
        sys.exit("Error! The following files did NOT have their domains predicted: " + ", ".join(domtableFiles - domtableBases))
    
    # find already processed files (assuming that if the pfd file exists, the pfs should too)
    alreadyDone = set()
    if not options.force_hmmscan:
        for domtable in domtableFiles:
            outputbase = ".".join(domtable.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(pfd_folder, outputbase + '.pfd')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(domtable)
    domtableFilesUnprocessed = domtableFiles - alreadyDone
    if len(domtableFilesUnprocessed) == 0: # Re-run
        print(" All domtable files had already been processed")
    elif len(alreadyDone) > 0: # Incomplete run
        if len(domtableFilesUnprocessed) < 20:
            print(" Warning! The following domtable files had not been processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in domtableFilesUnprocessed)))
        else:
            print(" Warning: {} domtable files will be processed".format(str(len(domtableFilesUnprocessed))))
    else: # First run
        print(" Processing {} domtable files".format(str(len(domtableFiles))))

    # If using the multiprocessing version and outputbase doesn't have any
    #  predicted domains, it's not as easy to remove if from the analysis
    #  (probably because parseHmmScan only has a copy of clusters et al?)
    # Using serialized version for now. Probably doesn't have too bad an impact
    #pool = Pool(cores,maxtasksperchild=32)
    for domtableFile in domtableFiles - alreadyDone:
        parseHmmScan(domtableFile, pfd_folder, pfs_folder, options.domain_overlap_cutoff)
        #pool.apply_async(parseHmmScan, args=(domtableFile,output_folder,options.domain_overlap_cutoff))
    #pool.close()
    #pool.join()
    
    # If number of pfd files did not change, no new sequences were added to the 
    #  domain fastas and we could try to resume the multiple alignment phase
    # baseNames have been pruned of BGCs with no domains that might've been added temporarily
    try_MA_resume = False
    if len(baseNames - set(pfd.split(os.sep)[-1][:-9] for pfd in alreadyDone)) == 0:
        try_MA_resume = True
    else:
        # new sequences will be added to the domain fasta files. Clean domains folder
        # We could try to make it so it's not necessary to re-calculate every alignment,
        #  either by expanding previous alignment files or at the very least, 
        #  re-aligning only the domain files of the newly added BGCs
        print(" New domain sequences to be added; cleaning domains folder")
        for thing in os.listdir(domains_folder):
            os.remove(os.path.join(domains_folder,thing))

    print(" Finished generating generating pfs and pfd files.")
    

    ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per sample objects
    print("\nProcessing domains sequence files")
    
    # All available pfd files
    allPfdFiles = set(glob(os.path.join(pfd_folder,"*.pfd")))
    
    # pfdFiles: all pfd files corresponding to the input files
    # (some input files could've been removed due to not having predicted domains)
    pfdFiles = set()
    for name in baseNames:
        pfdFiles.add(os.path.join(pfd_folder, name+".pfd"))
    
    # pfdBases: the actual set of input files that have pfd files
    pfdBases = allPfdFiles.intersection(pfdFiles)
    
    # verify previous step. 
    # All BGCs without predicted domains should no longer be in baseNames    
    if len(pfdFiles - pfdBases) > 0:
        sys.exit("Error! The following files did NOT have their domtable files processed: " + ", ".join(pfdFiles - pfdBases))

    filtered_matrix = []
    if options.skip_ma:
        print(" Running with skip_ma parameter: Assuming that the domains folder has all the fasta files")
        try:
            with open(os.path.join(cache_folder, "BGCs.dict"), "r") as BGC_file:
                BGCs = pickle.load(BGC_file)
                BGC_file.close()
        except IOError:
            sys.exit("BGCs file not found...")
    else:
        print(" Adding sequences to corresponding domains file")
            
        for outputbase in baseNames:
            if verbose:
                print("   Processing: " + outputbase)

            pfdFile = os.path.join(pfd_folder, outputbase + ".pfd")
            filtered_matrix = [[part.strip() for part in line.split('\t')] for line in open(pfdFile)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(bgc_fasta_folder, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and 
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not try_MA_resume:
                with open(fasta_file, "r") as fasta_file_handle:
                    fasta_dict = fasta_parser(fasta_file_handle) # all fasta info from a BGC
                save_domain_seqs(filtered_matrix, fasta_dict, domains_folder, outputbase)

            BGCs[outputbase] = BGC_dic_gen(filtered_matrix)
            
            del filtered_matrix[:]
            
        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(cache_folder, "BGCs.dict"), "wb") as BGC_file:
            pickle.dump(BGCs, BGC_file)
            BGC_file.close()
            
    # if it's a re-run, the pfd/pfs files were not changed, so the skip_ma flag
    # is activated. We have to open the pfd files to get the gene labels for
    # each domain
    # We now always have to have this data so the alignments are produced
    for outputbase in baseNames:
        DomainCountGene[outputbase] = array('B')
        corebiosynthetic_position[outputbase] = array('B')
        BGCGeneOrientation[outputbase] = array('b')
        pfdFile = os.path.join(pfd_folder, outputbase + ".pfd")
        filtered_matrix = [[part.strip() for part in line.split('\t')] for line in open(pfdFile)]
        
        domain_counter = 0
        gene_number = 0
        gene_label = filtered_matrix[0][-1] # initialize with first label
        has_corebio = False
        for row in filtered_matrix:
            if row[-1] != gene_label:
                # we changed to a new gene. Check whether previous has a 
                # core biosynthetic / anchor domain hit
                if has_corebio:
                    corebiosynthetic_position[outputbase].append(gene_number)
                    has_corebio = False
                
                if gene_label[-1] == "+":
                    BGCGeneOrientation[outputbase].append(1)
                else:
                    BGCGeneOrientation[outputbase].append(-1)
                    
                gene_label = row[-1] # update current label
                gene_number += 1 # advance gene number
                
                # record number of domains in previous gene
                DomainCountGene[outputbase].append(domain_counter)
                domain_counter = 1 # reset domain counter
            else:
                domain_counter += 1 # increase domain counter
                
            # TODO: if len(corebiosynthetic_position[outputbase]) == 0
            # do something with the list of pfam ids. Specifically, mark
            # (in this case TODO or always?) as biosynthetic genes, the ones that contain
            # domains from a special list. This list of special domains
            # comes from predicted domains within the CDSs marked as 'sec_met'
            # by antismash
            if row[-1] in bgc_info[outputbase].biosynthetic_genes:
                has_corebio = True
            
        # There is no transition when we finish, so analyze last gene
        if gene_label[-1] == "+":
            BGCGeneOrientation[outputbase].append(1)
        else:
            BGCGeneOrientation[outputbase].append(-1)
        DomainCountGene[outputbase].append(domain_counter)
        if has_corebio:
            corebiosynthetic_position[outputbase].append(gene_number)

    # Get the ordered list of domains
    print(" Reading the ordered list of domains from the pfs files")
    for outputbase in baseNames:
        pfsfile = os.path.join(pfs_folder, outputbase + ".pfs")
        if os.path.isfile(pfsfile):
            DomainList[outputbase] = get_domain_list(pfsfile)
        else:
            sys.exit(" Error: could not open " + outputbase + ".pfs")
                
                
    ### Step 5: Create SVG figures
    print(" Creating arrower-like figures for each BGC")
    
    # read hmm file. We'll need that info anyway for final visualization
    print("  Parsing hmm file for domain information")
    pfam_info = {}
    with open(os.path.join(pfam_dir, "Pfam-A.hmm"), "r") as pfam:
        putindict = False
        # assuming that the order of the information never changes
        for line in pfam:
            if line[:4] == "NAME":
                name = line.strip()[6:]
            if line[:3] == "ACC":
                acc = line.strip()[6:].split(".")[0]
            if line[:4] == "DESC":
                desc = line.strip()[6:]
                putindict = True
                
            if putindict:
                putindict = False
                pfam_info[acc] = (name, desc)
    print("    Done")
    
    # verify if there are figures already generated
    
    # All available SVG files
    availableSVGs = set()
    for svg in glob(os.path.join(svg_folder,"*.svg")):
        (root, ext) = os.path.splitext(svg)
        availableSVGs.add(root.split(os.sep)[-1])
        
    # Which files actually need to be generated
    working_set = baseNames - availableSVGs
    
    if len(working_set) > 0:
        color_genes = read_color_genes_file()
        color_domains = read_color_domains_file()
        pfam_domain_categories = read_pfam_domain_categories()
        
        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        print("  Reading BGC information and writing SVG")
        for bgc in working_set:
            SVG(False, os.path.join(svg_folder,bgc+".svg"), genbankDict[bgc][0], os.path.join(pfd_folder,bgc+".pfd"), True, color_genes, color_domains, pfam_domain_categories, pfam_info, bgc_info[bgc].records, bgc_info[bgc].max_width)
            
        color_genes.clear()
        color_domains.clear()
        pfam_domain_categories.clear()
    elif len(working_set) == 0:
        print("  All SVG from the input files seem to be in the SVG folder")
    
    availableSVGs.clear()
    print(" Finished creating figures")
    
    
    print("\n\n   - - Calculating distance matrix - -")
   
    # Do multiple alignments if needed
    if not options.skip_ma:
        print("Performing multiple alignment of domain sequences")
        
        # obtain all fasta files with domain sequences
        domain_sequence_list = set(glob(os.path.join(domains_folder,"*.fasta")))
        
        # compare with .algn set of files. Maybe resuming is possible if
        # no new sequences were added
        if try_MA_resume:
            temp_aligned = set(glob(os.path.join(domains_folder, "*.algn")))
            
            if len(temp_aligned) > 0:
                print(" Found domain fasta files without corresponding alignments")
                
                for a in temp_aligned:
                    if os.path.getsize(a) > 0:
                        domain_sequence_list.remove(a[:-5]+".fasta")
            
            temp_aligned.clear()
        
        # Try to further reduce the set of domain fastas that need alignment
        sequence_tag_list = set()
        header_list = []
        domain_sequence_list_temp = domain_sequence_list.copy()
        for domain_file in domain_sequence_list_temp:
            domain_name = ".".join(domain_file.split(os.sep)[-1].split(".")[:-1])
            
            # fill fasta_dict...
            with open(domain_file, "r") as fasta_handle:
                header_list = get_fasta_keys(fasta_handle)
                
            # Get the BGC name from the sequence tag. The form of the tag is:
            # >BGCXXXXXXX_BGCXXXXXXX_ORF25:gid...
            sequence_tag_list = set(s.split("_ORF")[0] for s in header_list)

            # ...to find out how many sequences do we actually have
            if len(sequence_tag_list) == 1:
                # avoid multiple alignment if the domains all belong to the same BGC
                domain_sequence_list.remove(domain_file)
                if verbose:
                    print(" Skipping Multiple Alignment for {} (appears only in one BGC)".format(domain_name))
        
        sequence_tag_list.clear()
        del header_list[:]
        
        domain_sequence_list_temp.clear()
            
        # Do the multiple alignment
        stop_flag = False
        if len(domain_sequence_list) > 0:
            print("\n Using hmmalign")
            launch_hmmalign(cores, domain_sequence_list)
                
            # verify all tasks were completed by checking existance of alignment files
            for domain_file in domain_sequence_list:
                if not os.path.isfile(domain_file[:-6]+".algn"):
                    print("   ERROR, {}.algn could not be found (possible issue with aligner).".format(domain_file[:-6]))
                    stop_flag = True
            if stop_flag:
                sys.exit()
                       
        else:
            print(" No domain fasta files found to align")
    
    
    # If there's something to analyze, load the aligned sequences
    if options_samples or options_all:
        print(" Trying to read domain alignments (*.algn files)")
        aligned_files_list = glob(os.path.join(domains_folder, "*.algn"))
        if len(aligned_files_list) == 0:
            sys.exit("No aligned sequences found in the domain folder (run without the --skip_ma parameter or point to the correct output folder)")
        for aligned_file in aligned_files_list:
            with open(aligned_file, "r") as aligned_file_handle:
                fasta_dict = fasta_parser(aligned_file_handle)
                for header in fasta_dict:
                    AlignedDomainSequences[header] = fasta_dict[header]

    clusterNames = tuple(sorted(clusters))

    # create output directory for network files
    network_files_folder = os.path.join(network_folder, run_name)
    create_directory(network_files_folder, "Network Files", False)

    # copy html templates
    dir_util.copy_tree(os.path.join(os.path.realpath(os.path.dirname(__file__)), "html_template", "output"), output_folder)

    # make a new run folder in the html output & copy the overview_html
    network_html_folder = os.path.join(output_folder, "html_content", "networks", run_name)
    create_directory(network_html_folder, "Network HTML Files", False)
    shutil.copy(os.path.join(os.path.realpath(os.path.dirname(__file__)), "html_template", "overview_html"), os.path.join(network_html_folder, "overview.html"))

    # create pfams.js
    pfams_js_file = os.path.join(output_folder, "html_content", "js", "pfams.js")
    if not os.path.isfile(pfams_js_file):
        with open(pfams_js_file, "w") as pfams_js:
            pfam_json = {}
            pfam_colors = generatePfamColorsMatrix(os.path.join(os.path.dirname(os.path.realpath(__file__)), "domains_color_file.tsv"))
            for pfam_code in pfam_info:
                pfam_obj = {}
                if pfam_code in pfam_colors:
                    pfam_obj["col"] = pfam_colors[pfam_code]
                else:
                    pfam_obj["col"] = "255,255,255"
                pfam_obj["desc"] = pfam_info[pfam_code][1]
                pfam_json[pfam_code] = pfam_obj
            pfams_js.write("var pfams={};\n".format(json.dumps(pfam_json, indent=4, separators=(',', ':'), sort_keys=True)))

    # track the sub modules (i.e. allOthers, allNRPS, ...)
    html_subs = []


    # Try to make default analysis using all files found inside the input folder
    if options_all:
        print("\nGenerating distance network files with ALL available input files")
    
        # This version contains *all* bgcs, including the ones from excluded classes
        # TODO: BiG-SCAPE class should contain ALL possible classes the bgc might
        # end up in
        print("   Writing the complete Annotations file for the complete set")
        network_annotation_path = os.path.join(network_files_folder, "Network_Annotations_Full.tsv")
        with open(network_annotation_path, "w") as network_annotation_file:
            network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
            for bgc in clusterNames:
                product = bgc_info[bgc].product
                network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
    
        # Making network files mixing all classes
        if options_mix:
            print("\n Mixing all BGC classes")
            
            # only choose from valid classes
            mix_set = []
            
            # Indexes ALL cluster names in the working set
            for clusterIdx,clusterName in enumerate(clusterNames):
                product = bgc_info[clusterName].product
                predicted_class = sort_bgc(product)
                if predicted_class.lower() in valid_classes:
                    mix_set.append(clusterIdx)
            
            # create output directory   
            create_directory(os.path.join(network_files_folder, "mix"), "  Mix", False)
            
            # Create an additional file with the list of all clusters in the class + other info
            print("   Writing annotation file")
            network_annotation_path = os.path.join(network_files_folder, "mix", "Network_Annotations_mix.tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in mix_set:
                    bgc = clusterNames[idx]
                    product = bgc_info[bgc].product
                    network_annotation_file.write("\t".join([bgc, 
                        bgc_info[bgc].accession_id, bgc_info[bgc].description, 
                        product, sort_bgc(product), bgc_info[bgc].organism, 
                        bgc_info[bgc].taxonomy]) + "\n")
                           
                
            print("  Calculating all pairwise distances")
            pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
            
            cluster_pairs = [(x, y, -1) for (x, y) in pairs]
            pairs.clear()
            network_matrix_mix = generate_network(cluster_pairs, cores)
            del cluster_pairs[:]
                
            print("  Writing output files")
            pathBase = os.path.join(network_files_folder, "mix")
            filenames = []
            for cutoff in cutoff_list:
                filenames.append(os.path.join(pathBase, "mix_c{:.2f}.network".format(cutoff)))
            cutoffs_and_filenames = list(zip(cutoff_list, filenames))
            del filenames[:]
            write_network_matrix(network_matrix_mix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
               
            print("  Calling Gene Cluster Families")
            reduced_network = []
            pos_alignments = {}
            for row in network_matrix_mix:
                reduced_network.append([int(row[0]), int(row[1]), row[2]])
                reverse = False
                if row[-1] == 1.0:
                    reverse = True
                pa = pos_alignments.setdefault(int(row[0]),{})
                # lcsStartA, lcsStartB, seedLength, reverse={True,False}
                pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
            del network_matrix_mix[:]
            html_subs.append({ "name" : "mix", "css" : "Others", "label" : "Mixed"})
            clusterJsonBatch(mix_set, pathBase, "mix", reduced_network, pos_alignments,
                             cutoffs=cutoff_list, clusterClans=options.clans,
                             clanCutoff=options.clan_cutoff, htmlFolder=network_html_folder)
            del mix_set[:]
            del reduced_network[:]
            
            
        # Making network files separating by BGC class
        if options_classify:
            print("\n Working for each BGC class")
            
            # reinitialize BGC_classes to make sure the bgc lists are empty
            BGC_classes = defaultdict(list)
        
            # Preparing gene cluster classes
            print("  Sorting the input BGCs\n")
            
            # Indexes ALL cluster names in the working set
            for clusterIdx,clusterName in enumerate(clusterNames):
                product = bgc_info[clusterName].product
                predicted_class = sort_bgc(product)
                if predicted_class.lower() in valid_classes:
                    BGC_classes[predicted_class].append(clusterIdx)
                
                # possibly add hybrids to 'pure' classes
                if options.hybrids:
                    if predicted_class == "PKS-NRP_Hybrids":
                        if "nrps" in valid_classes:
                            BGC_classes["NRPS"].append(clusterIdx)
                        if "t1pks" in product and "pksi" in valid_classes:
                            BGC_classes["PKSI"].append(clusterIdx)
                        if "t1pks" not in product and "pksother" in valid_classes:
                            BGC_classes["PKSother"].append(clusterIdx)
                    
                    if predicted_class == "Others" and "-" in product:
                        subclasses = set()
                        for subproduct in product.split("-"):
                            subclass = sort_bgc(subproduct)
                            if subclass.lower() in valid_classes:
                                subclasses.add(subclass)
                                
                        # Prevent mixed BGCs with sub-Others annotations to get
                        # added twice (e.g. indole-cf_fatty_acid has already gone
                        # to Others at this point)
                        if "Others" in subclasses:
                            subclasses.remove("Others")
                            
                            
                        for subclass in subclasses:
                            BGC_classes[subclass].append(clusterIdx)
                        subclasses.clear()

            # only make folders for the BGC_classes that are found
            for bgc_class in BGC_classes:
                print("\n  {} ({} BGCs)".format(bgc_class, str(len(BGC_classes[bgc_class]))))
                
                # create output directory   
                create_directory(os.path.join(network_files_folder, bgc_class), "  All - " + bgc_class, False)
                
                # Create an additional file with the final list of all clusters in the class
                print("   Writing annotation files")
                network_annotation_path = os.path.join(network_files_folder, bgc_class, "Network_Annotations_" + bgc_class + ".tsv")
                with open(network_annotation_path, "w") as network_annotation_file:
                    network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                    for idx in BGC_classes[bgc_class]:
                        bgc = clusterNames[idx]
                        product = bgc_info[bgc].product
                        network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
                    
                if len(BGC_classes[bgc_class]) > 1:
                    print("   Calculating all pairwise distances")
                    pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                    cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                    pairs.clear()
                    network_matrix = generate_network(cluster_pairs, cores)
                    del cluster_pairs[:]
                        
                    print("   Writing output files")
                    pathBase = os.path.join(network_files_folder, bgc_class)
                    filenames = []
                    for cutoff in cutoff_list:
                        filenames.append(os.path.join(pathBase, "{}_c{:.2f}.network".format(bgc_class, cutoff)))
                    cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                    del filenames[:]
                    write_network_matrix(network_matrix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)

                    print("  Calling Gene Cluster Families")
                    reduced_network = []
                    pos_alignments = {}
                    for row in network_matrix:
                        reduced_network.append([int(row[0]), int(row[1]), row[2]])
                        reverse = False
                        if row[-1] == 1.0:
                            reverse = True
                        pa = pos_alignments.setdefault(int(row[0]),{})
                        # lcsStartA, lcsStartB, seedLength, reverse={True,False}
                        pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
                    del network_matrix[:]

                    html_subs.append({ "name" : bgc_class, "css" : bgc_class, "label" : bgc_class})
                    clusterJsonBatch(BGC_classes[bgc_class], pathBase, bgc_class,
                                     reduced_network, pos_alignments, cutoffs=cutoff_list, 
                                     clusterClans=options.clans, clanCutoff=options.clan_cutoff, 
                                     htmlFolder=network_html_folder)
                    del BGC_classes[bgc_class][:]
                    del reduced_network[:]

    ## Try to make analysis for each sample
    if options_samples:
        pass
        #network_matrix_sample = []
        #clusterNames2idx = dict(zip(clusterNames,range(len(clusterNames))))
        
        #if len(sampleDict) == 1 and options_all:
            #print("\nNOT generating networks per sample (only one sample, covered in the all-vs-all case)")
        #else:
            #print("\nGenerating distance network files for each sample")

            #for sample, sampleClusters in sampleDict.items():
                #print("\n Sample: " + sample)
                #if len(sampleClusters) == 1:
                    #print(" Warning: Sample size = 1 detected. Not generating network for this sample ({})".format(sample))
                #else:
                    ## create output directory for this sample
                    #create_directory(os.path.join(network_files_folder, sample), " Samples - " + sample, False)
                        
                    ## Making network files mixing all classes
                    #if options_mix:
                        #print("\n  Mixing all BGC classes")
                        
                        #mix_set = []
                        #for clusterIdx, sampleCluster in enumerate(sampleClusters):
                            #product = bgc_info[sampleCluster].product
                            #predicted_class = sort_bgc(product)
                            #if predicted_class.lower() in valid_classes:
                                #mix_set.append(clusterIdx)
                        
                        ## Create an additional file with the list of all clusters in the class + other info
                        #print("   Writing annotation files")
                        #path_list = os.path.join(network_files_folder, sample, "Network_Annotations_Sample_" + sample + "_mix.tsv")
                        #with open(path_list, "w") as list_file:
                            #list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                            #for idx in mix_set:
                                #bgc = clusterNames[idx]
                                #product = bgc_info[bgc].product
                                #list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
            
                        #pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
                        
                        #cluster_pairs = [(x, y, -1) for (x, y) in pairs]
                        #pairs.clear()
                        #network_matrix_sample = generate_network(cluster_pairs, cores)
                        #del cluster_pairs[:]

                        #print("   Writing output files")
                        #pathBase = os.path.join(network_files_folder, sample, "sample_" + sample + "_mix")
                        #filenames = []
                        #for cutoff in cutoff_list:
                            #filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
                        #cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                        #write_network_matrix(network_matrix_sample, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
                        
                        #print("  Calling Gene Cluster Families")
                        #reduced_network = []
                        #for row in network_matrix_sample:
                            #reduced_network.append([int(row[0]), int(row[1]), row[2]])
                        #del network_matrix_sample[:]

                        #html_subs.append({ "name" : "sample_" + sample + "_mix", "css" : "Others", "label" : "S-" + sample + "-Mixed"})
                        #clusterJsonBatch(mix_set, pathBase, reduced_network, subName="sample_" + sample + "_mix", cutoffs=cutoff_list,clusterClans=options.clans,clanCutoff=options.clan_cutoff,htmlFolder=network_html_folder)
                        #del mix_set[:]
                        #del reduced_network[:]
                        
                    ## Making network files separating by BGC class
                    #if options_classify:
                        #print("\n  Working for each BGC class")

                        ## reinitialize BGC_classes to make sure the bgc lists are empty
                        #BGC_classes = defaultdict(list)
                    
                        ## Preparing gene cluster classes
                        #print("   Sorting the input BGCs\n")
                        #for cluster in sampleClusters:
                            #product = bgc_info[cluster].product
                            #predicted_class = sort_bgc(product)
                            #if predicted_class.lower() in valid_classes:
                                #BGC_classes[predicted_class].append(clusterNames2idx[cluster])
                            
                            ## possibly add hybrids to 'pure' classes
                            #if options.hybrids:
                                #if predicted_class == "PKS-NRP_Hybrids":
                                    #if "nrps" in valid_classes:
                                        #BGC_classes["NRPS"].append(clusterNames2idx[cluster])
                                    #if "t1pks" in product and "pksi" in valid_classes:
                                        #BGC_classes["PKSI"].append(clusterNames2idx[cluster])
                                    #if "t1pks" not in product and "pksother" in valid_classes:
                                        #BGC_classes["PKSother"].append(clusterNames2idx[cluster])
                                
                                #if predicted_class == "Others" and "-" in product:
                                    #subclasses = set()
                                    #for subproduct in product.split("-"):
                                        #subclass = sort_bgc(subproduct)
                                        #if subclass.lower() in valid_classes:
                                            #subclasses.add(subclass)
                                        
                                    ## Prevent mixed BGCs with sub-Others annotations to get
                                    ## added twice (e.g. indole-cf_fatty_acid has already gone
                                    ## to Others at this point)
                                    #if "Others" in subclasses:
                                        #subclasses.remove("Others")
                                        
                                    #for subclass in subclasses:
                                        #BGC_classes[subclass].append(clusterNames2idx[cluster])
                                    #subclasses.clear()

                        #for bgc_class in BGC_classes:
                            #folder_name = bgc_class
                                
                            #print("\n   " + folder_name + " (" + str(len(BGC_classes[bgc_class])) + " BGCs)")
                            #network_matrix_sample = []
                            
                            ## create output directory
                            #create_directory(os.path.join(network_files_folder, sample, folder_name), "   Sample " + sample + " - " + bgc_class, False)

                            ## Create an additional file with the final list of all clusters in the class
                            #print("   Writing annotation files")
                            #path_list = os.path.join(network_files_folder, sample, folder_name, "Network_Annotations_Sample_" + sample + "_" + folder_name + ".tsv")
                            #with open(path_list, "w") as list_file:
                                #list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                                #for idx in BGC_classes[bgc_class]:
                                    #bgc = clusterNames[idx]
                                    #product = bgc_info[bgc].product
                                    #list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")

                            #if len(BGC_classes[bgc_class]) > 1:
                                #pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                                
                                #cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                                #pairs.clear()
                                #network_matrix_sample = generate_network(cluster_pairs, cores)
                                #del cluster_pairs[:]
                                #print("    Writing output files")
                                #pathBase = os.path.join(network_files_folder, sample, folder_name,
                                                        #"sample_" + sample + "_" + folder_name)
                                #filenames = []
                                #for cutoff in cutoff_list:
                                    #filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
                                #cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                                #write_network_matrix(network_matrix_sample, cutoffs_and_filenames, include_singletons, clusterNames,bgc_info)
                                
                                #print("  Calling Gene Cluster Families")
                                #reduced_network = []
                                #for row in network_matrix_sample:
                                    #reduced_network.append([int(row[0]), int(row[1]), row[2]])
                                #del network_matrix_sample[:]
                                
                                #html_subs.append({ "name" : "sample_" + sample + "_" + folder_name, "css" : folder_name, "label" : "S-" + sample + "-" + folder_name})
                                #clusterJsonBatch(BGC_classes[bgc_class], pathBase, reduced_network, subName="sample_" + sample + "_" + folder_name, cutoffs=cutoff_list,clusterClans=options.clans,clanCutoff=options.clan_cutoff,htmlFolder=network_html_folder,className=bgc_class)
                                
                                #del BGC_classes[bgc_class][:]
                                #del reduced_network[:]

    # update bgc_results.js
    add_to_bigscape_results_js(run_name, html_subs, os.path.join(output_folder, "html_content", "js", "bigscape_results.js"))

    pickle.dump(bgc_info,open(os.path.join(cache_folder,'bgc_info.dict'),'wb'))
    runtime = time.time()-time1
    runtime_string = "\n\n\tMain function took {:.3f} s".format(runtime)
    with open(os.path.join(log_folder, "runtimes.txt"), 'a') as timings_file:
        timings_file.write(runtime_string + "\n")
    print(runtime_string)
    
