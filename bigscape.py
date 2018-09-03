#!/usr/bin/env python


"""
BiG-SCAPE

PI: Marnix Medema               marnix.medema@wur.nl

Main developers:
Jorge Navarro                   j.navarro@westerdijkinstitute.nl
Emmanuel (Emzo) de los Santos   E.De-Los-Santos@warwick.ac.uk
Satria Kautsar                  satria.kautsar@wur.nl


Usage:   Please see `python bigscape.py -h`

Example: python bigscape.py -c 8 --pfam_dir ./ -i ./inputfiles -o ./results

Status: beta

Official repository:
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
from itertools import product as combinations_product
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from argparse import ArgumentParser
from difflib import SequenceMatcher
from operator import itemgetter
import zipfile

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
import json
import shutil
from distutils import dir_util
from sklearn.cluster import AffinityPropagation
import networkx as nx


global use_relevant_mibig
global mibig_set
global genbankDict
global valid_classes

def process_gbk_files(gbk, min_bgc_size, bgc_info, files_no_proteins, files_no_biosynthetic_genes):
    """ Given a file path to a GenBank file, reads information about the BGC"""

    biosynthetic_genes = set()
    product_list_per_record = []
    fasta_data = []
    save_fasta = True
    adding_sequence = False
    contig_edge = False
    total_seq_length = 0
    record_end = 0
    offset_record_position = 0
    bgc_locus_tags = []
    locus_sequences = {}
    locus_coordinates = {}
    
    file_folder, fname = os.path.split(gbk)
    clusterName = fname[:-4]

    # See if we need to keep the sequence
    # (Currently) we have to open the file anyway to read all its 
    # properties for bgc_info anyway...
    outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
    if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0 and not force_hmmscan:
        if verbose:
            print(" File {} already processed".format(outputfile))
        save_fasta = False
    else:
        save_fasta = True
    
    try:
        # basic file verification. Substitutes check_data_integrity
        records = list(SeqIO.parse(gbk, "genbank"))
    except ValueError as e:
        print("   Error with file {}: \n    '{}'".format(gbk, str(e)))
        print("    (This file will be excluded from the analysis)")
        return
    else:
        total_seq_length = 0
        bgc_size = 0
        cds_ctr = 0
        product = "no type"
        del product_list_per_record[:]
        offset_record_position = 0
        
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
                        gene_id = CDS.qualifiers.get('gene',"")[0]
                        
                    
                    protein_id = ""
                    if "protein_id" in CDS.qualifiers:
                        protein_id = CDS.qualifiers.get('protein_id',"")[0]
                    
                    # nofuzzy_start/nofuzzy_end are obsolete
                    # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#nofuzzy_start
                    gene_start = offset_record_position + max(0, int(CDS.location.start))
                    gene_end = offset_record_position + max(0, int(CDS.location.end))
                    record_end = gene_end
                    
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
                                print("Warning, CDS ({}, {}) has fuzzy\
                                    start and end positions, and a \
                                    sequence length not multiple of \
                                    three. Skipping".format(clusterName, 
                                    CDS.qualifiers.get('locus_tag',"")[0]))
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
                
                
                    bgc_locus_tags.append(fasta_header)
                    locus_sequences[fasta_header] = prot_seq
                    locus_coordinates[fasta_header] = (gene_start, gene_end, len(prot_seq))
                    

            # TODO: if len(biosynthetic_genes) == 0, traverse record again
            # and add CDS with genes that contain domains labeled sec_met
            # we'll probably have to have a list of domains if we allow
            # fasta files as input
            
            # make absolute positions for ORFs in next records
            offset_record_position += record_end + 100
        
        if bgc_size > min_bgc_size:  # exclude the bgc if it's too small
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
                
            # Don't keep this bgc if its type not in valid classes specified by user
            # This will avoid redundant tasks like domain detection
            subproduct = set()
            for p in product.split("-"):
                subproduct.add(sort_bgc(p).lower())
            if "nrps" in subproduct and ("pksi" in subproduct or "pksother" in subproduct):
                subproduct.add("pks-nrp_hybrids")
                
            if len(valid_classes & subproduct) == 0:
                if verbose:
                    print(" Skipping {} (type: {})".format(clusterName, product))
                return False
            
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
            bgc_info[clusterName] = bgc_data(records[0].id, records[0].description, product, len(records), max_width, bgc_size + (record_count-1)*1000, records[0].annotations["organism"], ",".join(records[0].annotations["taxonomy"]), biosynthetic_genes.copy(), contig_edge)

            if len(bgc_info[clusterName].biosynthetic_genes) == 0:
                files_no_biosynthetic_genes.append(clusterName+".gbk")

            # TODO why re-process everything if it was already in the list?
            # if name already in genbankDict.keys -> add file_folder
            # else: extract all info
            if clusterName in genbankDict.keys():
                # Name was already in use. Use file_folder as the new sample's name
                genbankDict[clusterName][1].add(file_folder) 
            else:
                # See if we need to write down the sequence
                if total_seq_length > 0:
                    # location of first instance of the file is genbankDict[clustername][0]
                    genbankDict.setdefault(clusterName, [gbk, set([file_folder])])

                    if save_fasta:
                        # Find overlaps in CDS regions and delete the shortest ones.
                        # This is thought as a solution for selecting genes with 
                        # alternate splicing events
                        # Food for thought: imagine CDS A overlapping CDS B overlapping
                        # CDS C. If len(A) > len(B) > len(C) and we first compare A vs B
                        # and delete A, then B vs C and delete B: would that be a better
                        # solution than removing B? Could this actually happen?
                        # TODO What if the overlapping CDS is in the reverse strand?
                        #  maybe it should be kept as it is
                        # TODO what are the characterized differences in prokarytote
                        #  vs eukaryote CDS overlap?
                        del_list = set()
                        for a, b in combinations(bgc_locus_tags, 2):
                            a_start, a_end, a_len = locus_coordinates[a]
                            b_start, b_end, b_len = locus_coordinates[b]
                            
                            if b_end <= a_start or b_start >= a_end:
                                pass
                            else:
                                # calculate overlap
                                if a_start > b_start:
                                    ov_start = a_start
                                else:
                                    ov_start = b_start

                                if a_end < b_end:
                                    ov_end = a_end
                                else:
                                    ov_end = b_end

                                overlap_length = ov_end - ov_start
                                
                                # allow the overlap to be as large as 10% of the
                                # shortest CDS. Overlap length is in nucleotides
                                # here, whereas a_len, b_len are protein 
                                # sequence lengths
                                if overlap_length/3 > 0.1*min(a_len,b_len):
                                    if a_len > b_len:
                                        del_list.add(b)
                                    else:
                                        del_list.add(a)
                        
                        for locus in del_list:
                            if verbose:
                                print("   Removing {} because it overlaps with other ORF".format(locus))
                            bgc_locus_tags.remove(locus)
                        
                        with open(outputfile,'w') as fastaHandle:
                            for locus in bgc_locus_tags:
                                fastaHandle.write("{}\n".format(locus))
                                fastaHandle.write("{}\n".format(locus_sequences[locus]))
                            adding_sequence = True
                else:
                    files_no_proteins.append(fname)

            if verbose:
                print("  Adding {} ({} bps)".format(fname, str(bgc_size)))
                                
        else:
            print(" Discarding {} (size less than {} bp, was {})".format(clusterName, str(min_bgc_size), str(bgc_size)))
    
    return adding_sequence


def get_gbk_files(inputpath, outputdir, bgc_fasta_folder, min_bgc_size, exclude_gbk_str, bgc_info):
    """Searches given directory for genbank files recursively, will assume that
    the genbank files that have the same name are the same genbank file. 
    Returns a dictionary that contains the names of the clusters found as keys
    and a list that contains [0] a path to the genbank file and [1] the 
    samples that the genbank file is a part of.
    Extract and write the sequences as fasta files if not already in the Fasta 
    folder.
    return: {cluster_name:[genbank_path,[s_a,s_b...]]}
    """
    file_counter = 0
    processed_sequences = 0
    files_no_proteins = []
    files_no_biosynthetic_genes = []


    if os.path.isfile(inputpath):
        files = [inputpath]
    else:
        # Unfortunately, this does not work in Python 2:
        #files = glob(os.path.join(inputpath,"**/*.gbk"), recursive=True) 
        files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(inputpath)
                 for f in files if f.endswith(".gbk")]
        
    for filepath in files:
        file_folder, fname = os.path.split(filepath)
        
        if type(exclude_gbk_str) == str and exclude_gbk_str != "" and \
                                            exclude_gbk_str in fname:
            if verbose:
                print(" Skipping file " + fname)
            continue
        elif type(exclude_gbk_str) == list and exclude_gbk_str != [] and \
                        any([word in fname for word in exclude_gbk_str]):
            print(" Skipping file " + fname)
            continue
        if "_ORF" in fname:
            print(" Skipping file {} (string '_ORF' is used internally)".format(fname))
            continue
        
        if " " in filepath:
            sys.exit("\nError: Input GenBank files should not have spaces in their path as hmmscan cannot process them properly ('too many arguments').")
        
        file_counter += 1
        if process_gbk_files(filepath, min_bgc_size, bgc_info, files_no_proteins, files_no_biosynthetic_genes):
            processed_sequences += 1
    
    if len(files_no_proteins) > 0:
        print("  Warning: Input set has files without protein sequences. They will be discarded")
        print("   (See no_sequences_list.txt)")
        with open(os.path.join(outputdir, "no_sequences_list.txt"), "w") as noseqs:
            for f in sorted(files_no_proteins):
                noseqs.write("{}\n".format(f))
        
    if len(files_no_biosynthetic_genes) > 0 and (mode == "glocal" or mode == "auto"):
        print("  Warning: Input set has files with no Biosynthetic Genes (affects alignment mode)")
        print("   See no_biosynthetic_genes_list.txt")
        with open(os.path.join(outputdir, "logs", "no_biosynthetic_genes_list.txt"), "w") as nobiogenes:
            for f in sorted(files_no_biosynthetic_genes):
                nobiogenes.write("{}\n".format(f))
    
    print("\n Starting with {:d} files".format(file_counter))
    print(" Files that had its sequence extracted: {:d}".format(processed_sequences))

    return


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
    
    dist, jaccard, dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB, seedLength, reverse = cluster_distance_lcs(cluster1, cluster2, domain_list_A,
        domain_list_B, dcg_a, dcg_b, core_pos_a, core_pos_b, go_a, go_b, bgc_class)
        
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

        return 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, S, S_anchor, 0, 0, 0, 0


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
    if s > sr or (s == sr and go_A[a] == go_b[b]):
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
        
    elif s < sr:
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
        
    # special cases: s == sr
    
    # if only one gene matches, choose the one with the most domains
    elif s == 1:
        seqmatch = SequenceMatcher(None, A_string, b_string)
        max_domains = 0
        x = 0   # index in A with the gene with most domains
        y = 0   # index in B with the gene with most domains
        for a, b, z in seqmatch.get_matching_blocks():
            if z != 0:
                if DomainCountGene[A][a] > max_domains:
                    x = a
                    y = b
                    max_domains = DomainCountGene[A][a]
        
        # note: these slices are in terms of genes, not domains (which are 
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        sliceStartA = x
        sliceLengthA = 1
        sliceLengthB = 1
        
        if BGCGeneOrientation[A][x] == BGCGeneOrientation[B][y]:
            sliceStartB = y
            dcg_B = dcg_b
            B_string = b_string
            reverse = False
            b_name = B
        else:
            sliceStartB = len(BGCGeneOrientation[B]) - y - 1
            dcg_B = list(reversed(dcg_b))
            B_string = b_string_reverse
            reverse = True
            b_name = B + "*"
               
    # s == sr and (s == 0 or s > 1)
    else:
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

    if mode=="glocal" or (mode=="auto" and (bgc_info[A].contig_edge or bgc_info[B].contig_edge)):
        #X: bgc that drive expansion
        #Y: the other bgc
        # forward: True if expansion is to the right
        # returns max_score, final positions for X and Y
        #score_expansion(X_string, dcg_X, Y_string, dcg_Y, downstream=True/False)
            
        # Expansion is relatively costly. We ask for a minimum of 3 genes
        # for the core overlap before proceeding with expansion.
        biosynthetic_hit_A = False
        for biosynthetic_position in core_pos_A:
            if biosynthetic_position >= sliceStartA and biosynthetic_position <= (sliceStartA+sliceLengthA):
                biosynthetic_hit_A = True
                break
        if sliceLengthA >= 3 or biosynthetic_hit_A:
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
    
    domain_file_stk = domain_file[:-6]+".stk"
    hmmalign_pars = ["hmmalign", "-o", domain_file_stk, "-", domain_file]
    proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout, stdout=subprocess.PIPE, shell=False)
    
    proc_hmmfetch.stdout.close()
    proc_hmmalign.communicate()[0]
    proc_hmmfetch.wait()
    
    if verbose:
        print(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))
    
    stockholm_parser(domain_file_stk)
    #SeqIO.convert(domain_file_stk, "stockholm", domain_file[:-6]+".algn", "fasta")
    
    
def stockholm_parser(stkFile):
    reference = ""
    algnDict = defaultdict(str)
    algnFile = stkFile[:-3] + 'algn'
    if not os.path.isfile(algnFile): # prevents overwriting algn files
        with open(stkFile, 'r') as infile:
            for l in infile:
                line = l.strip()
                if line.startswith("#=GC RF"):
                    reference += line[7:].strip()
                elif line == "":
                    continue
                elif line[0] == "/" or line[0] == "#":
                    continue
                else:
                    a = line.split(" ")
                    header = a[0]
                    algn = a[-1]
                    algnDict[header] += algn
                    
        # get start-end coordinates of every "x" island (original consensus)
        # in the reference
        state_reference = False
        slicing_tuples = []
        for pos in range(len(reference)):
            if reference[pos] == "x" and not state_reference:
                state_reference = True
                start = pos
            if reference[pos] == "." and state_reference:
                state_reference = False
                slicing_tuples.append((start,pos))
        if state_reference:
            slicing_tuples.append((start, len(reference)))
    
    if len(algnDict) > 0:
        with open(algnFile, "w") as outfile:
            for header in algnDict:
                sequence = ""
                for a,b in slicing_tuples:
                    sequence += algnDict[header][a:b]
                outfile.write(">{}\n".format(header))
                outfile.write(sequence + "\n")
    return


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
            if outputbase in mibig_set:
                mibig_set.remove(outputbase)
            
    else:
        sys.exit("Error: hmmscan file " + outputbase + " was not found! (parseHmmScan)")

    return("")


def clusterJsonBatch(bgcs, pathBase, className, matrix, pos_alignments, cutoffs=[1.0], damping=0.9, clusterClans=False, clanCutoff=(0.5,0.8), htmlFolder=None):
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

    family_data = { # will be returned, for use in overview.js
        "label": className,
        "families": [],
        "families_similarity": []
    }
    
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
    # If needed, bgcInt2Ext[i] would simply be bgcs[i]
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
        bgcJsonDict[bgcName]["end"] = bgc_info[bgcName].bgc_size
        bgcJsonDict[bgcName]["mibig"] = bgcName in mibig_set
        
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
        bs_data_js.write("dataLoaded('bs_data');\n")
    shutil.copy(os.path.join(os.path.realpath(os.path.dirname(__file__)), "html_template", "index_html"), os.path.join(module_html_path, "index.html"))
    
    # Create network
    g = nx.Graph()
    
    for cutoff in cutoffs:
        print("  Cutoff: {}".format(cutoff))
        # first task is to find labels (exemplars) for each node.
        # Assign disconnected (singleton) BGCs to themselves
        labels = [0 for i in range(numBGCs)]
        
        # clear all edges from network
        g.clear()
        
        # add all nodes
        g.add_nodes_from(bgcs)
        
        # add all (allowed) edges
        for bgc1 in bgcs:
            for bgc2 in simDict.get(bgc1,{}).keys():
                if simDict[bgc1][bgc2] > 1 - cutoff:
                    g.add_edge(bgc1, bgc2)
                    
        for subgraph in nx.connected_components(g):
            numBGCs_subgraph = len(subgraph)
            # smaller subgraphs don't need to be clustered
            if numBGCs_subgraph < 3:
                temp = list(subgraph)
                for bgc in temp:
                    labels[bgcExt2Int[bgc]] = bgcExt2Int[temp[0]]
            else:
                bgcExt2Sub_ = dict(zip([c for c in subgraph], range(numBGCs_subgraph)))
                bgcSub2Ext_ = dict(zip(range(numBGCs_subgraph), [c for c in subgraph]))
                                   
                simMatrix = np.zeros((numBGCs_subgraph, numBGCs_subgraph), dtype=np.float32)
                for bgc1, bgc2 in combinations(subgraph,2):
                    try:
                        simMatrix[bgcExt2Sub_[bgc1], bgcExt2Sub_[bgc2]] = simDict[bgc1][bgc2]
                    except KeyError:
                        simMatrix[bgcExt2Sub_[bgc1], bgcExt2Sub_[bgc2]] = simDict[bgc2][bgc1]
                        
                af = AffinityPropagation(damping=damping, max_iter=1000, convergence_iter=200, affinity="precomputed").fit(simMatrix)
                labelsSub = af.labels_
                exemplarsSub = af.cluster_centers_indices_
                
                # TODO go straight to familiesDict
                for i in range(numBGCs_subgraph):
                    labels[bgcExt2Int[bgcSub2Ext_[i]]] = bgcExt2Int[bgcSub2Ext_[exemplarsSub[labelsSub[i]]]]
        
        # Recalculate distance matrix as we'll need it with clans
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
            print("   ...done")

        bs_distances = [[float("{:.3f}".format(simMatrix[row, col])) for col in 
                         range(row+1)] for row in range(numBGCs)]
        
        
        familiesDict = defaultdict(list)
        for i in range(numBGCs):
            familiesDict[bgcs[labels[i]]].append(bgcs[i])
        familyIdx = sorted(familiesDict.keys()) # identifiers for each family
        

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
                newick_trees[exemplar_idx] = "({}):0.01;".format(",".join([str(bgcExt2Int[x])+":0.0" for x in gcf]))
                
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
                if verbose:
                    print("  Warning: core shared domains for GCF {} consists of a single domain ({})".format(exemplar_idx, [x for x in tree_domains][0]))
            
            # Get the alignments of the core domains
            alignments = {}
            # initialize every sequence alignment entry. Don't do defaultdict!
            alignments[exemplar_idx] = ""
            
            out_of_tree_bgcs = [] # bgcs that don't share a common domain core
            delete_list = set() # remove this bgcs from alignment
            gcf.remove(exemplar_idx) # separate exemplar from the rest of the bgcs
            for bgc in gcf:
                alignments[bgc] = ""

            match_dict = {}
            for domain in tree_domains:
                specific_domain_list_A = BGCs[exemplar][domain]
                num_copies_a = len(specific_domain_list_A)
                for exemplar_domain_copy in specific_domain_list_A:
                    alignments[exemplar_idx] += AlignedDomainSequences[exemplar_domain_copy]
                
                seq_length = len(AlignedDomainSequences[specific_domain_list_A[0]])
                
                for bgc in alignments:
                    match_dict.clear()
                    if bgc == exemplar_idx:
                        pass
                    elif domain not in domain_sets[bgc]:
                        if verbose:
                            print("   BGC {} ({}) does not share a common domain core (GCF: {}, domain: {})".format(clusterNames[bgc], bgc, exemplar_idx, domain))
                        out_of_tree_bgcs.append(bgc)
                        delete_list.add(bgc)
                    elif bgc in delete_list:
                        pass
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
                                
                                matches = 0
                                gaps = 0
                                
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
                    
                for bgc in list(delete_list):
                    del alignments[bgc]
                delete_list.clear()
                
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
                
            # make tree
            newick_file_path = os.path.join(gcf_trees_path, "GCF_c{:4.2f}_{:05d}.newick".format(cutoff,exemplar_idx))
            with open(newick_file_path, "w") as newick_file:
                command = ["fasttree", "-nopr", "-quiet", alignment_file_path]
                p = subprocess.Popen(command, stdout=newick_file, shell=False)
                p.wait() # only with process has terminated will the file be ready

            # read tree, post-process it and save it
            if not os.path.isfile(newick_file_path) or os.path.getsize(newick_file_path) == 0:
                print(newick_file_path)
                sys.exit(" ERROR: newick file not created or empty (GCF_c{:4.2f}_{:05d})".format(cutoff,exemplar_idx))
            else:
                with open(newick_file_path,"r") as newick_file:
                    try:
                        tree = Phylo.read(newick_file, 'newick')
                    except ValueError as e:
                        print(" Warning! There was an error while reading tree file {}".format(newick_file))
                        print(str(e))
                        newick_trees[exemplar_idx] = ""
                    else:
                        try:
                            tree.root_at_midpoint()
                        except UnboundLocalError:
                            # Noticed this could happen if the sequences are exactly
                            # the same and all distances == 0
                            if verbose:
                                print(" Warning: Unable to root at midpoint file {}".format(newick_file_path))
                            pass
                        newick = tree.format("newick")
                        
                        # convert branches' names to indices for visualization
                        for name in bgc_name_to_idx:
                            newick = newick.replace(name, str(bgcExt2Int[bgc_name_to_idx[name]]))

                        newick_trees[exemplar_idx] = newick
       
        ### - - - GCC - - -
        bs_similarity_families = []
        if clusterClans and cutoff == clanClassificationCutoff:
            # Detect if there's only 1 GCF. It makes pySAPC crash
            if len(familyIdx) == 1:
                clanLabels = [1]
                continue
            
            #famSimMatrix = lil_matrix((len(familyIdx), len(familyIdx)), dtype=np.float32)
            famSimMatrix = np.zeros((len(familyIdx), len(familyIdx)), dtype=np.float32)
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
                    
            # if we have the identity matrix here, it means all GCFs are separate
            # (nothing to cluster). Note: still can crash if values are 
            # sufficiently low. Catch this error later
            #if np.count_nonzero(simMatrix) == 0:
                #clanLabels = []
                #continue
            
            # add main diagonal
            for family in range(len(familyIdx)):
                famSimMatrix[family,family] = 1.0
                
            bs_similarity_families = famSimMatrix.tolist()
            
            #clanLabels = pysapc.SAP(damping=damping, max_iter=500,
                                #preference='min').fit_predict(famSimMatrix)
            af = AffinityPropagation(damping=damping, max_iter=1000, convergence_iter=200, affinity="precomputed").fit(famSimMatrix)
            labelsClans = af.labels_
            exemplarsClans = af.cluster_centers_indices_
            
            # affinity propagation can fail in some circumstances (e.g. only singletons)
            if exemplarsClans is not None:
                # translate and record GCF label instead of GCF number
                clanLabels = [familyIdx[exemplarsClans[labelsClans[i]]] for i in range(len(familyIdx))]
            else:
                clanLabels = []
            
        else:
            clanLabels = []

        if len(clanLabels) > 0:
            clansDict = defaultdict(list)
            for i in range(len(familyIdx)):
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
        
        family_data["families"] = []
        for family, members in familiesDict.items():
            family_data["families"].append({
                "label": "FAM_{:05d}".format(family),
                "members": members # use external indexing
            })

        # Positional alignment information is based on DomainCountGene, which
        # does not contain empty genes (i.e. with no domains). 
        domainGenes2allGenes = {}
        
        ## BGC Family alignment information
        bs_families_alignment = []
        for family, members in familiesDict.items():
            for bgc in members:
                domainGenes2allGenes[bgc] = {}
                has_domains = 0
                for orf in range(len(bs_data[bgcExt2Int[bgc]]["orfs"])):
                    if len(bs_data[bgcExt2Int[bgc]]["orfs"][orf]["domains"]) > 0:
                        domainGenes2allGenes[bgc][has_domains] = orf
                        has_domains += 1
                        
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
                        
                        if length == 0:
                            pass
                        elif reverse:
                            # special case. bgc was reference (first) in lcs
                            a = domainGenes2allGenes[family][len(DomainCountGene[clusterNames[family]])-a-length]
                            b = domainGenes2allGenes[bgc][b+length-1] # -1 go to 0-index
                        else:
                            a = domainGenes2allGenes[family][a]
                            b = domainGenes2allGenes[bgc][b]
                    else:
                        a = domainGenes2allGenes[family][a]
                        if length == 0:
                            pass
                        
                        elif reverse:
                            
                            b = domainGenes2allGenes[bgc][len(DomainCountGene[clusterNames[bgc]])-b-1]
                        else:
                            b = domainGenes2allGenes[bgc][b]
                    
                    
                    if length == 0:
                        length = 1
                        # let's try aligning using the genes with most domains
                        # after all, they ended up being in the same GCF
                        # for some reason
                        x = max(DomainCountGene[clusterNames[family]])
                        x = DomainCountGene[clusterNames[family]].index(x)
                        a = domainGenes2allGenes[family][x]
                        
                        y = max(list(DomainCountGene[clusterNames[bgc]]))
                        y = DomainCountGene[clusterNames[bgc]].index(y)
                        
                        #check orientation
                        if BGCGeneOrientation[clusterNames[family]][x] == BGCGeneOrientation[clusterNames[bgc]][y]:
                            b = domainGenes2allGenes[bgc][y]
                            reverse = False
                        else:
                            b = domainGenes2allGenes[bgc][len(DomainCountGene[clusterNames[bgc]])-y-1]
                            reverse = True
                            
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

        # get family-family similarity matrix
        bs_similarity_families = [[get_composite_bgc_similarities([bgcs[bid] for bid in bs_families[row]["members"]], [bgcs[bid] for bid in bs_families[col]["members"]], simDict) if (row != col) else (1.00, (1.00, bgcs[bs_families[row]["members"][0]], bgcs[bs_families[row]["members"][0]]), (1.00, bgcs[bs_families[row]["members"][0]], bgcs[bs_families[row]["members"][0]])) for col in 
                         range(row+1)] for row in range(len(bs_families))]

        family_data["families_similarity"] = bs_similarity_families;

        ## Write bgc_networks.js
        with open(os.path.join(module_html_path, "bs_networks.js"), "w") as bs_networks_js:
            bs_networks_js.write("var bs_similarity={};\n".format(json.dumps(bs_distances, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families={};\n".format(json.dumps(bs_families, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families_alignment={};\n".format(json.dumps(bs_families_alignment, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_similarity_families={};\n".format(json.dumps(bs_similarity_families, indent=4, separators=(',', ':'), sort_keys=True)))
            if len(clanLabels) > 0:
                bs_networks_js.write("var bs_clans={};\n".format(json.dumps(bs_clans, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("dataLoaded('bs_networks');\n")

        
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
                    
    return family_data


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
    
    parser.add_argument("-l", "--label", dest="label", help="An extra label for\
                        this run (will be used as part of the folder name within\
                        the network_files results)")
    
    parser.add_argument("-i", "--inputdir", dest="inputdir", 
                        default=os.path.dirname(os.path.realpath(__file__)),
                        help="Input directory of gbk files, if left empty, all \
                        gbk files in current and lower directories will be used.")

    parser.add_argument("--exclude_gbk_str", dest="exclude_gbk_str", 
                        default="final", nargs="+",
                        help="If any string in this list occurs in the gbk \
                        filename, this file will not be used for the analysis.\
                        (default: final)")
    
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

    parser.add_argument("--cutoffs", dest="cutoffs", nargs="+", default=[0.30], 
                        type=float, choices=[FloatRange(0.0, 1.0)], 
                        help="Generate networks using multiple raw distance \
                    cutoff values, example: --cutoffs 0.1, 0.25, 0.5, 1.0. Default: \
                    c=0.3.")
                    
    parser.add_argument("--clans-off", dest="clans",action="store_false", 
                        default=True, help="BiG-SCAPE will perform a second \
                        layer of clustering and attempt to group families \
                        assigned from clustering with cutoff of 0.5 to clans")

    parser.add_argument("--clan_cutoff",dest="clan_cutoff",default=[0.3,0.7], 
                        type=float, choices=[FloatRange(0.0, 1.0)],nargs=2,
                        help="Cutoff Parameters for which clustering families \
                        into clans will be performed in raw distance. First \
                        value is the cutoff value family assignments for BGCs \
                        used in clan clustering (default: 0.3). Second value is \
                        the cutoff value for clustering families into clans \
                        (default: 0.7). Average linkage for BGCs in a family is\
                        used for distances between families. Example: \
                        --clan_cutoff 0.3 0.7)")

    parser.add_argument("--hybrids-off", dest="hybrids", action="store_false", 
                        default=True, help="Toggle to also add BGCs with hybrid\
                        predicted products from the PKS/NRPS Hybrids and Others\
                        classes to each subclass (e.g. a 'terpene-nrps' BGC from\
                        Others would be added to the Terpene and NRPS classes)")
    
    parser.add_argument("--mode", dest="mode", default="glocal", choices=["global",
                            "glocal", "auto"], help="Alignment mode for each pair of\
                            gene clusters. 'global': the whole list of domains \
                            of each BGC are compared; 'glocal': Longest Common \
                            Subcluster mode. Redefine the subset of the domains \
                            used to calculate distance by trying to find the \
                            longest slice of common domain content per gene in \
                            both BGCs, then expand each slice. 'auto': use glocal\
                            when at least one of the BGCs in each pair has the \
                            'contig_edge' annotation from antiSMASH v4+, \
                            otherwise use global mode on that pair")
    
    parser.add_argument("--anchorfile", dest="anchorfile", 
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),"anchor_domains.txt"),
                        help="Provide a custom location for the anchor domains \
                        file, default is anchor_domains.txt.")
                    
    parser.add_argument("--force_hmmscan", dest="force_hmmscan", action="store_true", 
                        default=False, help="Force domain prediction using \
                        hmmscan even if BiG-SCAPE finds processed domtable files\
                        (e.g. to use a new version of PFAM).")
    parser.add_argument("--skip_ma", dest="skip_ma", action="store_true", 
                        default=False, help="Skip multiple alignment of domains'\
                        sequences. Use if alignments have been generated in a \
                        previous run.")
    
    parser.add_argument("--mibig", dest="use_relevant_mibig", action=
        "store_true", default=False, help="Use included BGCs from then MIBiG \
        database. Only relevant (i.e. those with distance < max(cutoffs) against\
        the input set) will be used. Using version (version 1.3). See https://mibig.secondarymetabolites.org/")
     
    parser.add_argument("--query_bgc", help="Instead of making an all-VS-all \
                    comparison of all the input BGCs, choose one BGC to \
                    compare with the rest of the set (one-VS-all). The \
                    query BGC does not have to be within inputdir")
    
    parser.add_argument("--domain_whitelist", help="Only analyze include those\
                        BGCs that include domains with the pfam accessions \
                        found in the domain_whitelist.txt file", default=False,
                        action="store_true")

    parser.add_argument("--version", action="version", version="%(prog)s 201809")

    return parser.parse_args()


if __name__=="__main__":
    options = CMD_parser()
    
    class bgc_data:
        def __init__(self, accession_id, description, product, records, max_width, bgc_size, organism, taxonomy, biosynthetic_genes, contig_edge):
            # These two properties come from the genbank file:
            self.accession_id = accession_id
            self.description = description
            # AntiSMASH predicted class of compound:
            self.product = product
            # number of records in the genbank file (think of multi-locus BGCs):
            self.records = records
            # length of largest record (it will be used for ArrowerSVG):
            self.max_width = int(max_width)
            # length of the entire bgc (can include several records/subclusters)
            self.bgc_size = bgc_size
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
    
    global force_hmmscan
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
    global run_data
    run_data = {}

    global clusterNames, bgcClassNames
    
    include_singletons = options.include_singletons
    
    cores = int(options.cores)
    
    options_mix = options.mix
    options_classify = not options.no_classify
    
    force_hmmscan = options.force_hmmscan
    
    mode = options.mode
    
    cutoff_list = options.cutoffs
    for c in cutoff_list:
        if c <= 0.0 or c > 1.0:
            print(" Removing invalid cutoff value {}".format(str(c)))
            cutoff_list.remove(c)
    max_cutoff = max(cutoff_list)
            
    # if we want to classify by clans make sure that the clanCutoff is included
    # in the cutoffs to do AP in clusterJsonBatch
    if options.clans:
        fc, cc = options.clan_cutoff
        if fc not in cutoff_list:
            if fc <= 0.0 or fc > 1.0:
                sys.exit("Error: invalid cutoff value for GCF calling")
            else:
                cutoff_list = sorted(cutoff_list.append(fc))
            
        if cc <= 0.0 or cc > 1.0:
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

    has_query_bgc = False
    if options.query_bgc:
        has_query_bgc = True
        if not os.path.isfile(options.query_bgc):
            sys.exit("Error: Query BGC not found")
            
    verbose = options.verbose
    
    run_mode_string = ""
    networks_folder_all = "networks_all"
    if options.hybrids:
        networks_folder_all += "_hybrids"
        run_mode_string += "_hybrids"
    if mode == "auto":
        networks_folder_all += "_auto"
        run_mode_string += "_auto"
    elif mode == "glocal":
        networks_folder_all += "_glocal"
        run_mode_string += "_glocal"
    else:
        run_mode_string += "_global"
    
    time1 = time.time()

    start_time = time.localtime()
    run_name = "{}{}".format(time.strftime("%Y-%m-%d_%H-%M-%S", start_time), run_mode_string)
    if options.label:
        run_name = run_name + "_" + options.label
    run_data["start_time"] = time.strftime("%d/%m/%Y %H:%M:%S", start_time)
    run_data["parameters"] = " ".join(sys.argv[1:])
    run_data["input"] = {}

    # Make the following available for possibly deleting entries within parseHmmScan
    global gbk_files, sampleDict, clusters, baseNames
    
    
    # Get domain_whitelist
    has_whitelist = False
    if options.domain_whitelist:
        bigscape_path = os.path.dirname(os.path.realpath(__file__))
        if os.path.isfile(os.path.join(bigscape_path,"domain_whitelist.txt")):
            domain_whitelist = set()
            for line in open(os.path.join(bigscape_path,"domain_whitelist.txt"), "r"):
                if line[0] == "#":
                    continue
                domain_whitelist.add(line.split("\t")[0].strip())
            if len(domain_whitelist) == 0:
                print("Error: --domain_whitelist used, but no domains found in the file")
            else:
                has_whitelist = True
        else:
            sys.exit("Error: domain_whitelist.txt file not found")
    
    
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
    
    # finally, define weights for mix
    bgc_class_weight["mix"] = (0.2, 0.75, 0.05, 2.0) # default when not separating in classes
    BGC_classes = defaultdict(list)
    # mix class will always be the last element of the tuple
    bgcClassNames = tuple(sorted(list(bgc_class_weight)) + ["mix"])
    assert bgcClassNames[-1] == 'mix'


    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    bgc_info = {} # Stores, per BGC: predicted type, gbk Description, number of records, width of longest record, GenBank's accession, Biosynthetic Genes' ids
    genbankDict = {}
    
    # Exclude single string
    exclude_gbk_str = options.exclude_gbk_str
    if type(exclude_gbk_str) == str and exclude_gbk_str != "":
        print(" Skipping files with '{}' in their filename".format(exclude_gbk_str))
    # Exclude from list of strings
    elif type(exclude_gbk_str) == list and exclude_gbk_str != []:
        print(" Skipping files with one or more of the following strings in \
            their filename: {}".format(", ".join(exclude_gbk_str)))
    
    # Read included MIBiG
    # Change this for every officially curated MIBiG bundle
    # (file, final folder, number of bgcs)
    mibig_zipfile_numbgcs = ("MIBiG_1.3_gbks.zip", "1.3+_final_gbks", 1393)
    use_relevant_mibig = options.use_relevant_mibig
    mibig_set = set()
    if use_relevant_mibig:
        print("\n Trying to read bundled MIBiG BGCs as reference")
        mibig_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"Annotated_MIBiG_reference")
        bgcs_path = os.path.join(mibig_path,mibig_zipfile_numbgcs[1])
        
        # try to see if the zip file has already been decompressed
        numbgcs = len(glob(os.path.join(bgcs_path,"*.gbk")))
        if numbgcs == 0:
            if not zipfile.is_zipfile(os.path.join(mibig_path,mibig_zipfile_numbgcs[0])):
                sys.exit("Did not find file {}. Please re-download it from the official repository".format(mibig_zipfile_numbgcs[0]))
                
            with zipfile.ZipFile(os.path.join(mibig_path,mibig_zipfile_numbgcs[0]), 'r') as mibig_zip:
                for fname in mibig_zip.namelist():
                    if fname[-3:] != "gbk":
                        continue
                
                    extractedbgc = mibig_zip.extract(fname,path=mibig_path)
                    if verbose:
                        print("  Extracted {}".format(extractedbgc))
        
        elif mibig_zipfile_numbgcs[2] == numbgcs:
            print("  MIBiG BGCs seem to have been extracted already")
        else:
            sys.exit("Did not find the correct number of MIBiG BGCs ({}). Please clean the 'Annotated MIBiG reference' folder from any .gbk files first".format(mibig_zipfile_numbgcs[2]))
        
        print("\nImporting MIBiG files")
        get_gbk_files(bgcs_path, output_folder, bgc_fasta_folder, int(options.min_bgc_size), exclude_gbk_str, bgc_info)
        
        for i in genbankDict.keys():
            mibig_set.add(i)
            
    
    print("\nImporting GenBank files")
    get_gbk_files(options.inputdir, output_folder, bgc_fasta_folder, int(options.min_bgc_size), exclude_gbk_str, bgc_info)
    
    if has_query_bgc:
        query_bgc = ".".join(options.query_bgc.split(os.sep)[-1].split(".")[:-1])
        if query_bgc in genbankDict:
            print("\nQuery BGC already added")
            pass
        else:
            print("\nImporting query BGC file")
            get_gbk_files(options.query_bgc, output_folder, bgc_fasta_folder, int(options.min_bgc_size), exclude_gbk_str, bgc_info)
            
        if query_bgc not in genbankDict:
            sys.exit("Error: not able to include Query BGC (check valid classes, BGC size, etc. Run again with --verbose)")
    # clusters and sampleDict contain the necessary structure for all-vs-all and sample analysis
    clusters = list(genbankDict.keys())
    
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
    if force_hmmscan:
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
                # verify domtable content
                with open(outputfile, "r") as domtablefile:
                    for line in domtablefile.readlines():
                        if line.startswith("# Option settings:"):
                            linecols = line.split()
                            if "hmmscan" in linecols and "--domtblout" in linecols:
                                alreadyDone.add(fasta)
                                break
                
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
    if not force_hmmscan:
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

    print(" Finished generating pfs and pfd files.")
    

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
    pfd_dict_domains = defaultdict(int)
    orf_keys = {}
    for outputbase in baseNames:
        DomainCountGene[outputbase] = array('B')
        corebiosynthetic_position[outputbase] = array('H')
        BGCGeneOrientation[outputbase] = array('b')
        pfdFile = os.path.join(pfd_folder, outputbase + ".pfd")
        
        #pfd_dict_domains contains the number of domains annotated in the
        # pfd file for each orf tag
        with open(pfdFile,"r") as pfdf:
            for line in pfdf:
                pfd_dict_domains[line.strip().split("\t")[-1]] += 1
        
        # extract the orf number from the tag and use it to traverse the BGC
        for orf in pfd_dict_domains.keys():
            orf_num = int(orf.split(":")[0].split("_ORF")[1])
            orf_keys[orf_num] = orf
        
        orf_num = 0
        for orf_key in sorted(orf_keys.keys()):
            orf = orf_keys[orf_key]
            if orf[-1] == "+":
                BGCGeneOrientation[outputbase].append(1)
            else:
                BGCGeneOrientation[outputbase].append(-1)
                
            DomainCountGene[outputbase].append(pfd_dict_domains[orf])
            
            if orf in bgc_info[outputbase].biosynthetic_genes:
                corebiosynthetic_position[outputbase].append(orf_num)
            orf_num += 1
        
        pfd_dict_domains.clear()
        orf_keys.clear()

        ## TODO: if len(corebiosynthetic_position[outputbase]) == 0
        ## do something with the list of pfam ids. Specifically, mark
        ## (in this case TODO or always?) as biosynthetic genes, the ones that contain
        ## domains from a special list. This list of special domains
        ## comes from predicted domains within the CDSs marked as 'sec_met'
        ## by antismash
            
            
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
        color_genes = {}
        color_domains = read_color_domains_file()
        pfam_domain_categories = {}
        
        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        print("  Reading BGC information and writing SVG")
        for bgc in working_set:
            with open(genbankDict[bgc][0],"r") as handle:
                SVG(False, os.path.join(svg_folder,bgc+".svg"), handle, bgc, os.path.join(pfd_folder,bgc+".pfd"), True, color_genes, color_domains, pfam_domain_categories, pfam_info, bgc_info[bgc].records, bgc_info[bgc].max_width)
        
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
    
    # we have to find the idx of query_bgc
    if has_query_bgc:
        try:
            query_bgc_idx = clusterNames.index(query_bgc)
        except ValueError:
            sys.exit("Error finding the index of Query BGC")

    # create output directory for network files
    network_files_folder = os.path.join(network_folder, run_name)
    create_directory(network_files_folder, "Network Files", False)
    run_data["networks"] = []

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
    print("\nGenerating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    print("   Writing the complete Annotations file for the complete set")
    network_annotation_path = os.path.join(network_files_folder, "Network_Annotations_Full.tsv")
    with open(network_annotation_path, "w") as network_annotation_file:
        network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
        for bgc in clusterNames:
            product = bgc_info[bgc].product
            network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
    
    
    # Find index of all MIBiG BGCs if necessary
    if use_relevant_mibig:
        name_to_idx = {}
        for clusterIdx,clusterName in enumerate(clusterNames):
            name_to_idx[clusterName] = clusterIdx
            
        mibig_set_indices = set()
        for bgc in mibig_set:
            mibig_set_indices.add(name_to_idx[bgc])

    # Making network files mixing all classes
    if options_mix:
        print("\n Mixing all BGC classes")
        
        # only choose from valid classes
        mix_set = []
        
        # create working set with indices of valid clusters
        for clusterIdx,clusterName in enumerate(clusterNames):
            if has_whitelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DomainList[clusterName]})
                    
                if len(domain_whitelist & bgc_domain_set) == 0:
                    continue
            
            product = bgc_info[clusterName].product
            predicted_class = sort_bgc(product)
            
            if predicted_class.lower() in valid_classes:
                mix_set.append(clusterIdx)
        
        print("\n  {} ({} BGCs)".format("Mix", str(len(mix_set))))

        # create output directory
        create_directory(os.path.join(network_files_folder, "mix"), "  Mix", False)
        
        print("  Calculating all pairwise distances")
        if has_query_bgc:
            pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx], mix_set)])
        else:
            # convert into a set of ordered tuples
            pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
        
        cluster_pairs = [(x, y, -1) for (x, y) in pairs]
        pairs.clear()
        network_matrix_mix = generate_network(cluster_pairs, cores)
        
        del cluster_pairs[:]

        # add the rest of the edges in the "Query network"
        if has_query_bgc:
            new_set = []
            
            # rows from the distance matrix that will be pruned 
            del_list = [] 
        
            for idx, row in enumerate(network_matrix_mix):
                a, b, distance = int(row[0]), int(row[1]), row[2]
                
                if a == b:
                    continue
                
                if distance <= max_cutoff:
                    if a == query_bgc_idx:
                        new_set.append(b)
                    else:
                        new_set.append(a)
                else:
                    del_list.append(idx)
            
            for idx in sorted(del_list, reverse=True):
                del network_matrix_mix[idx]
            del del_list[:]
            
            pairs = set([tuple(sorted(combo)) for combo in combinations(new_set, 2)])
            cluster_pairs = [(x, y, -1) for (x, y) in pairs]
            pairs.clear()
            network_matrix_new_set = generate_network(cluster_pairs, cores)
            del cluster_pairs[:]
            
            # Update the network matrix (QBGC-vs-all) with the distances of
            # QBGC's GCF
            network_matrix_mix.extend(network_matrix_new_set)
            
            # Update actual list of BGCs that we'll use
            mix_set = new_set
            mix_set.extend([query_bgc_idx])
            mix_set.sort() # clusterJsonBatch expects ordered indices
        
            # Create an additional file with the list of all clusters in the class + other info
            # This version of the file only has information on the BGCs connected to Query BGC
            print("   Writing annotation file")
            network_annotation_path = os.path.join(network_files_folder, "mix", "Network_Annotations_mix_QueryBGC.tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in mix_set:
                    bgc = clusterNames[idx]
                    product = bgc_info[bgc].product
                    network_annotation_file.write("\t".join([bgc, 
                        bgc_info[bgc].accession_id, bgc_info[bgc].description, 
                        product, sort_bgc(product), bgc_info[bgc].organism, 
                        bgc_info[bgc].taxonomy]) + "\n")
        elif use_relevant_mibig:
            n = nx.Graph()
            n.add_nodes_from(mix_set)
            mibig_set_del = []
            network_matrix_set_del = []
            
            for idx, row in enumerate(network_matrix_mix):
                a, b, distance = int(row[0]), int(row[1]), row[2]
                if distance <= max_cutoff:
                    n.add_edge(a, b, index=idx)
                    
            for component in nx.connected_components(n): # note: 'component' is a set
                numBGCs_subgraph = len(component)
                
                # catch if the subnetwork is comprised only of MIBiG BGCs
                if len(component & mibig_set_indices) == numBGCs_subgraph:
                    for bgc in component:
                        mibig_set_del.append(bgc)
                    
            # Get all edges between bgcs marked for deletion
            for (a, b, idx) in n.subgraph(mibig_set_del).edges.data('index'):
                network_matrix_set_del.append(idx)
                
            # delete all edges between marked bgcs
            for row_idx in sorted(network_matrix_set_del, reverse=True):
                del network_matrix_mix[row_idx]
            del network_matrix_set_del[:]
            
            print("   Removing {} non-relevant MIBiG BGCs".format(len(mibig_set_del)))
            mix_set_idx = 0
            bgc_to_mix_set_idx = {}
            for idx, bgc in enumerate(mix_set):
                bgc_to_mix_set_idx[bgc] = idx
            
            for bgc_idx in sorted(mibig_set_del, reverse=True):
                del mix_set[bgc_to_mix_set_idx[bgc_idx]]
            del mibig_set_del[:]
            

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
        family_data = clusterJsonBatch(mix_set, pathBase, "mix", reduced_network, pos_alignments,
                            cutoffs=cutoff_list, clusterClans=options.clans,
                            clanCutoff=options.clan_cutoff, htmlFolder=network_html_folder)
        run_data["networks"].append(family_data)
        del mix_set[:]
        del reduced_network[:]
        
        
    # Making network files separating by BGC class
    if options_classify:
        print("\n Working for each BGC class")
        
        # reinitialize BGC_classes to make sure the bgc lists are empty
        BGC_classes = defaultdict(list)
    
        # Preparing gene cluster classes
        print("  Sorting the input BGCs\n")
        
        # create and sort working set for each class
        for clusterIdx,clusterName in enumerate(clusterNames):
            if has_whitelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DomainList[clusterName]})
                    
                if len(domain_whitelist & bgc_domain_set) == 0:
                    continue
            
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
            if has_query_bgc:
                # not interested in this class if our Query BGC is not here...
                if query_bgc_idx not in BGC_classes[bgc_class]:
                    continue
            
            print("\n  {} ({} BGCs)".format(bgc_class, str(len(BGC_classes[bgc_class]))))
            if use_relevant_mibig:
                if len(set(BGC_classes[bgc_class]) & mibig_set_indices) == len(BGC_classes[bgc_class]):
                    print(" - All clusters in this class are MIBiG clusters -")
                    print("  If you'd like to analyze MIBiG clusters, turn off the --mibig option")
                    print("  and point --inputdir to the Annotated_MIBiG_reference folder")
                    continue
            
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
            
            print("   Calculating all pairwise distances")
            if has_query_bgc:
                pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx],BGC_classes[bgc_class])])
            else:
                pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                
            cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
            pairs.clear()
            network_matrix = generate_network(cluster_pairs, cores)
            #pickle.dump(network_matrix,open("others.ntwrk",'wb'))
            del cluster_pairs[:]
            #network_matrix = pickle.load(open("others.ntwrk", "rb"))
                
            # add the rest of the edges in the "Query network"
            if has_query_bgc:
                new_set = []
                
                # rows from the distance matrix that will be pruned 
                del_list = []
                
                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]
                    
                    # avoid QBGC-QBGC
                    if a == b:
                        continue
                    
                    if distance <= max_cutoff:
                        if a == query_bgc_idx:
                            new_set.append(b)
                        else:
                            new_set.append(a)
                    else:
                        del_list.append(idx)
                
                for idx in sorted(del_list, reverse=True):
                    del network_matrix[idx]
                del del_list[:]
                
                pairs = set([tuple(sorted(combo)) for combo in combinations(new_set, 2)])
                cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                pairs.clear()
                network_matrix_new_set = generate_network(cluster_pairs, cores)
                del cluster_pairs[:]
                                    
                # Update the network matrix (QBGC-vs-all) with the distances of
                # QBGC's GCF
                network_matrix.extend(network_matrix_new_set)
                
                # Update actual list of BGCs that we'll use
                BGC_classes[bgc_class] = new_set
                BGC_classes[bgc_class].extend([query_bgc_idx])
                BGC_classes[bgc_class].sort()
                
                # Create an additional file with the list of all clusters in the class + other info
                # This version of the file only has information on the BGCs connected to Query BGC
                print("   Writing annotation file (Query BGC)")
                network_annotation_path = os.path.join(network_files_folder, bgc_class, "Network_Annotations_" + bgc_class + "_QueryBGC.tsv")
                with open(network_annotation_path, "w") as network_annotation_file:
                    network_annotation_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                    for idx in BGC_classes[bgc_class]:
                        bgc = clusterNames[idx]
                        product = bgc_info[bgc].product
                        network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
            elif use_relevant_mibig:
                n = nx.Graph()
                n.add_nodes_from(BGC_classes[bgc_class])
                mibig_set_del = []
                network_matrix_set_del = []
                
                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]
                    if distance <= max_cutoff:
                        n.add_edge(a, b, index=idx)
                        
                for component in nx.connected_components(n): # note: 'component' is a set
                    numBGCs_subgraph = len(component)
                    
                    # catch if the subnetwork is comprised only of MIBiG BGCs
                    if len(component & mibig_set_indices) == numBGCs_subgraph:
                        for bgc in component:
                            mibig_set_del.append(bgc)
                
                # Get all edges between bgcs marked for deletion
                for (a, b, idx) in n.subgraph(mibig_set_del).edges.data('index'):
                    network_matrix_set_del.append(idx)
                
                # delete all edges between marked bgcs
                for row_idx in sorted(network_matrix_set_del, reverse=True):
                    del network_matrix[row_idx]
                del network_matrix_set_del[:]
                            
                print("   Removing {} non-relevant MIBiG BGCs".format(len(mibig_set_del)))
                bgc_to_class_idx = {}
                for idx, bgc in enumerate(BGC_classes[bgc_class]):
                    bgc_to_class_idx[bgc] = idx
                for bgc_idx in sorted(mibig_set_del, reverse=True):
                    del BGC_classes[bgc_class][bgc_to_class_idx[bgc_idx]]
                del mibig_set_del[:]
                    
                
                
            if len(BGC_classes[bgc_class]) < 2:
                continue
                
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
            family_data = clusterJsonBatch(BGC_classes[bgc_class], pathBase, bgc_class,
                                reduced_network, pos_alignments, cutoffs=cutoff_list, 
                                clusterClans=options.clans, clanCutoff=options.clan_cutoff, 
                                htmlFolder=network_html_folder)
            run_data["networks"].append(family_data)
            del BGC_classes[bgc_class][:]
            del reduced_network[:]

    # fetch genome list for overview.js
    genomes = []
    classes = []
    clusterNamesToGenomes = {}
    clusterNamesToClasses = {}
    inputClustersIdx = [] # contain only indexes (from clusterNames) of input BGCs (non-mibig)
    for idx, bgc in enumerate(clusterNames):
        if bgc in mibig_set:
            continue
        inputClustersIdx.append(idx)
        # get class info
        product = bgc_info[bgc].product
        predicted_class = sort_bgc(product)
        if predicted_class not in classes:
            clusterNamesToClasses[bgc] = len(classes)
            classes.append(predicted_class)
        else:
            clusterNamesToClasses[bgc] = classes.index(predicted_class)
        # get identifier info
        identifier = ""
        if len(bgc_info[bgc].organism) > 1:
            identifier = bgc_info[bgc].organism
        else : # use original genome file name (i.e. exclude "..clusterXXX from antiSMASH run")
            file_name_base = os.path.splitext(os.path.basename(genbankDict[bgc][0]))[0]
            identifier = file_name_base.rsplit(".cluster",1)[0]
        if len(identifier) < 1:
            identifier = "Unknown Genome {}".format(len(genomes))
        if identifier not in genomes:
            clusterNamesToGenomes[bgc] = len(genomes)
            genomes.append(identifier)
        else:
            clusterNamesToGenomes[bgc] = genomes.index(identifier)
    run_data["input"]["accession"] = [{ "id": "genome_{}".format(i), "label": acc } for i, acc in enumerate(genomes)]
    run_data["input"]["accession_newick"] = [] # todo ...
    run_data["input"]["classes"] = [{ "label": cl } for cl in classes ] # todo : colors
    run_data["input"]["bgc"] = [{ "id": clusterNames[idx], "acc": clusterNamesToGenomes[clusterNames[idx]], "class": clusterNamesToClasses[clusterNames[idx]] } for idx in inputClustersIdx]


    # update family data (convert global bgc indexes into input-only indexes)
    for network in run_data["networks"]:
        for family in network["families"]:
            new_members = []
            mibig = []
            for bgcIdx in family["members"]:
                if bgcIdx in inputClustersIdx:                    
                    new_members.append(inputClustersIdx.index(bgcIdx))
                else: # is a mibig bgc
                    clusterName = clusterNames[bgcIdx]
                    if clusterName in mibig_set:
                        mibig.append(clusterName)
            family["mibig"] = mibig
            family["members"] = new_members


    # generate overview data
    end_time = time.localtime()
    duration = int(time.mktime(end_time)) - int(time.mktime(start_time))
    run_data["end_time"] = time.strftime("%d/%m/%Y %H:%M:%S", end_time)
    run_data["duration"] = "{}h{}m{}s".format((duration // 3600), ((duration % 3600) // 60), ((duration % 3600) % 60))
    with open(os.path.join(network_html_folder, "run_data.js"), "w") as run_data_js:
        run_data_js.write("var run_data={};\n".format(json.dumps(run_data, indent=4, separators=(',', ':'), sort_keys=True)))
        run_data_js.write("dataLoaded();\n");

    # update bgc_results.js
    add_to_bigscape_results_js(run_name, html_subs, os.path.join(output_folder, "html_content", "js", "bigscape_results.js"))

    pickle.dump(bgc_info,open(os.path.join(cache_folder,'bgc_info.dict'),'wb'))
    runtime = time.time()-time1
    runtime_string = "\n\n\tMain function took {:.3f} s".format(runtime)
    with open(os.path.join(log_folder, "runtimes.txt"), 'a') as timings_file:
        timings_file.write(runtime_string + "\n")
    print(runtime_string)
    
