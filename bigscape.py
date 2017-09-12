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

from Bio import SeqIO
from Bio.SeqFeature import BeforePosition, AfterPosition
from Bio import AlignIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix

from functions import *
from ArrowerSVG import *

import numpy as np
from array import array
from scipy.sparse import lil_matrix
from scipy.optimize import linear_sum_assignment
import pysapc


def get_gbk_files(inputdir, bgc_fasta_folder, min_bgc_size, exclude_gbk_str, bgc_info):
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
    
    print("\nImporting GenBank files")
    if exclude_gbk_str != "":
        print(" Skipping files with '" + exclude_gbk_str + "' in their filename")

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
            
            if exclude_gbk_str != "" and exclude_gbk_str in fname:
                print(" Skipping file " + fname)
                continue
            if "_ORF" in fname:
                print(" Skipping file {} (string '_ORF' is used internally)".format(fname))
                continue
            
            if " " in fname:
                sys.exit("\nError: Input GenBank files should not have spaces in their filenames as HMMscan cannot process them properly ('too many arguments').")
            
            # See if we need to write down the sequence
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
                        if "cluster" in feature.type and "product" in feature.qualifiers:
                            if len(feature.qualifiers["product"]) > 1:
                                print("  WARNING: more than product annotated in record " + str(record_count) + ", " + fname)
                                break
                            else:
                                product_list_per_record.append(feature.qualifiers["product"][0].replace(" ",""))
                
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
                            

                            if save_fasta:
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
                                        
                                fasta_data.append((fasta_header, prot_seq))
                
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
                    bgc_info[clusterName] = bgc_data(records[0].id, records[0].description, product, len(records), max_width, biosynthetic_genes.copy(), contig_edge)

                    if len(bgc_info[clusterName].biosynthetic_genes) == 0:
                        print("  Warning: No Biosynthetic Genes found in " + clusterName)

                    # TODO why re-process everything if it was already in the list?
                    # if name already in genbankDict.keys -> add current_dir
                    # else: extract all info
                    if clusterName in genbankDict.keys():
                        # Name was already in use. Use current_dir as the new sample's name
                        genbankDict[clusterName][1].add(current_dir) 
                    else:
                        # location of first instance of the file is genbankDict[clustername][0]
                        genbankDict.setdefault(clusterName, [os.path.join(dirpath, fname), set([current_dir])])
                        
                        # See if we need to write down the sequence
                        outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
                        if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                            processed_sequences += 1
                        else:
                            with open(outputfile,'w') as fastaHandle:
                                for header_sequence in fasta_data:
                                    fastaHandle.write("{}\n".format(str(header_sequence[0])))
                                    fastaHandle.write("{}\n".format(str(header_sequence[1])))
                            
                            
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
    
    print("\n Starting with {} files".format(str(file_counter)))
    print(" Files that had its sequence extracted: " + str(file_counter - processed_sequences))

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
            with open(os.path.join(output_folder, "runtimes.txt"), 'a') as timings_file:
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

        # cluster1Idx, cluster2Idx, bgcClassIdx, distance, jaccard, DSS, AI, rDSSNa, rDSSa, S, Sa
        return array('f',[cluster1Idx,cluster2Idx,bgcClassIdx,  1,0,0,0,0,0,1,1])
    

    dist, jaccard, dss, ai, rDSSna, rDSS, S, Sa = cluster_distance(cluster1, cluster2,
                                                                   domain_list_A, domain_list_B, bgc_class) #sequence dist
        
    network_row = array('f',[cluster1Idx, cluster2Idx, bgcClassIdx, dist, (1-dist)**2, jaccard, dss, ai, rDSSna, rDSS, S, Sa])
    
    return network_row
    

def cluster_distance(a, b, a_domlist, b_domlist, bgc_class): 
    """Compare two clusters using information on their domains, and the 
    sequences of the domains"""

    Jaccardw, DSSw, AIw, anchorboost = bgc_class_weight[bgc_class]

    temp_domain_fastas = {}
    
    A = a
    B = b
    A_domlist = a_domlist[:]
    B_domlist = b_domlist[:]
    
    setA = set(A_domlist)
    setB = set(B_domlist)
    intersect = setA & setB
    
    S = 0
    S_anchor = 0
    
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
        
        return 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, S, S_anchor

    
    # define the subset of domain sequence tags to include in
    # the DSS calculation. This is done for every domain.
    A_domain_sequence_slice_bottom = defaultdict(int)
    A_domain_sequence_slice_top = defaultdict(int)
    B_domain_sequence_slice_bottom = defaultdict(int)
    B_domain_sequence_slice_top = defaultdict(int)
    
    
    # In local mode, try to align the shorter BGC ("BGC-fragment") to the
    # best matching slice of the larger BGC
    if local:
        # BGC A will be the shortest
        if len(a_domlist) < len(b_domlist):
            A = a
            B = b
            A_domlist = a_domlist[:]
            tmpB_domlist = b_domlist[:]
        else:
            A = b
            B = a
            A_domlist = b_domlist[:]
            tmpB_domlist = a_domlist[:]
        
        # Find the slice of the larger BGC where the shorter one fits the best
        setA = set(A_domlist)
        intersect = set()
        startB = 0
        lengthA = len(A_domlist) # length of raw list including copies
        for i in range(len(tmpB_domlist) - lengthA + 1):
            tmpBset = set(tmpB_domlist[i:i+lengthA])
            if len(setA.intersection(tmpBset)) > len(intersect):
                startB = i
                intersect = setA.intersection(tmpBset)
        B_domlist = tmpB_domlist[startB:startB+lengthA]
        setB = set(B_domlist)

        # initialize domain sequence slices
        for domain in setA:
            A_domain_sequence_slice_bottom[domain] = 0
            A_domain_sequence_slice_top[domain] = len(BGCs[A][domain])
            
        # the longest BGC needs to be sliced for domain copies as well
        for i in range(startB):
            domain = tmpB_domlist[i]
            B_domain_sequence_slice_bottom[domain] += 1
            
        # for each top, start at bottom
        for domain in setB:
            B_domain_sequence_slice_top[domain] = B_domain_sequence_slice_bottom[domain]
        for i in range(startB, startB+lengthA):
            domain = tmpB_domlist[i]
            B_domain_sequence_slice_top[domain] += 1

    else:
        # initialize domain sequence slices
        for domain in setA:
            A_domain_sequence_slice_bottom[domain] = 0
            A_domain_sequence_slice_top[domain] = len(BGCs[A][domain])
            
        for domain in setB:
            B_domain_sequence_slice_bottom[domain] = 0
            B_domain_sequence_slice_top[domain] = len(BGCs[B][domain])
        
    
    # JACCARD INDEX
    Jaccard = len(intersect)/ len(setA|setB)
    

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
            unshared_occurrences = BGCs[A][unshared_domain]
        except KeyError:
            unshared_occurrences = BGCs[B][unshared_domain]
            
        # don't look at domain version, hence the split
        if unshared_domain[:7] in anchor_domains:
            domain_difference_anchor += len(unshared_occurrences)
        else:
            domain_difference += len(unshared_occurrences)
                    
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
                            
                DistanceMatrix[domsa][domsb] = 1- ( matches/(seq_length-gaps) )
                
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
    
    if len(A_domlist) < 2 or len(B_domlist) < 2:
        AI = 0.0
    else:
        setA_pairs = set()
        for l in range(len(A_domlist)-1):
            setA_pairs.add(tuple(sorted([A_domlist[l],A_domlist[l+1]])))
        
        setB_pairs = set()
        for l in range(len(B_domlist)-1):
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
        
    return Distance, Jaccard, DSS, AI, DSS_non_anchor, DSS_anchor, S, S_anchor


"""
*FFT-NS-i (iterative refinement method; two cycles only):
mafft --retree 2 --maxiterate 2 input [> output]
fftnsi input [> output]
*FFT-NS-i (iterative refinement method; max. 1000 iterations):
mafft --retree 2 --maxiterate 1000 input [> output]
*FFT-NS-2 (fast; progressive method):
mafft --retree 2 --maxiterate 0 input [> output]
fftns input [> output]
*FFT-NS-1 (very fast; recommended for >2000 sequences; progressive method with a rough guide tree):
mafft --retree 1 --maxiterate 0 input [> output]
*NW-NS-i (iterative refinement method without FFT approximation; two cycles only):
mafft --retree 2 --maxiterate 2 --nofft input [> output]
nwnsi input [> output]
*NW-NS-2 (fast; progressive method without the FFT approximation):
mafft --retree 2 --maxiterate 0 --nofft input [> output]
nwns input [> output]
*NW-NS-PartTree-1 (recommended for ~10,000 to ~50,000 sequences; progressive method with the PartTree algorithm):
mafft --retree 1 --maxiterate 0 --nofft --parttree input [> output]

With FFT NS 1, a distance matrix is first generated using the 6 tuple score between each pair of sequences both sequences
are scanned from the start for matching 6 tuples, and when a match is found the score is incremented and scanning continues
from the next residue [4]. A guide tree is then constructed by clustering according to these distances, and the sequences 
are then aligned using the branching order of the guide tree. With FFT NS 2, the alignment produced by the FFT NS 1 method
is used to regenerate the distance matrix and the guide tree, and then do a second progressive alignment. In this paper,
FFT NS 1 will be specified whenever distance measures are needed. If no distance measures are required, the default 
FFTNS2 method will be used. 
"""
@timeit
def run_mafft(al_method, maxit, cores, mafft_pars, domain):
    """Runs mafft program with the provided parameters.
    The specific parameters used to run mafft with are actually very important for the final result.
    Using mafft with the most progressive parameters does indeed affect the quality of the distance.
    It is better to just use the domain information for the distance if more computationally intensive options
    for mafft cannot be used. Setting maxiterate to anything higher than 10 did not make a difference in accuracy in the nrps testset"""
    
    alignment_file = domain + ".algn"
    
    
    mafft_cmd_list = []
    mafft_cmd_list.append("mafft") 
    #mafft_cmd_list.append("--distout") #distout will save the distance matrix in a ".hat2" file
    mafft_cmd_list.append("--quiet")
    mafft_cmd_list.append(al_method)
    if maxit != 0:
        mafft_cmd_list.append("--maxiterate " + str(maxit))
        
    mafft_cmd_list.append("--thread")
    mafft_cmd_list.append(str(cores))
    
    if mafft_pars != "":
        mafft_cmd_list.append(mafft_pars)
        
    mafft_cmd_list.append(str(domain) + ".fasta")
    mafft_cmd_list.append(">")
    mafft_cmd_list.append(str(alignment_file))
    
    mafft_cmd = " ".join(mafft_cmd_list)
    
    if verbose:
        print("   " + mafft_cmd)
    subprocess.check_output(mafft_cmd, shell=True)


def launch_hmmalign(cores, domains):
    """
    Launches instances of hmmalign with multiprocessing.
    Note that the domains parameter contains the .fasta extension
    """
    pool = Pool(cores, maxtasksperchild=32)
    pool.map(run_hmmalign, domains)
    pool.close()
    pool.join()
    
def run_hmmalign(domain):
    #domain already contains the full path, with the file extension
    domain_base = domain.split(os.sep)[-1][:-6]
    hmmfetch_pars = ["hmmfetch", os.path.join(pfam_dir,"Pfam-A.hmm.h3m"), domain_base]
    proc_hmmfetch = subprocess.Popen(hmmfetch_pars, stdout=subprocess.PIPE, shell=False)
    
    hmmalign_pars = ["hmmalign", "-o", domain.replace(".fasta",".stk"), "-", domain]
    proc_hmmalign = subprocess.Popen(hmmalign_pars, stdin=proc_hmmfetch.stdout, stdout=subprocess.PIPE, shell=False)
    
    proc_hmmfetch.stdout.close()
    proc_hmmalign.communicate()[0]
    proc_hmmfetch.wait()
    
    if verbose:
        print(" ".join(hmmfetch_pars) + " | " + " ".join(hmmalign_pars))
    
    SeqIO.convert(domain[:-6]+".stk", "stockholm", domain[:-6]+".algn", "fasta")
    

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

def clusterJsonBatch(outputFileBase,matrix,cutoffs=[1.0],damping=0.8,clusterClans=False,clanCutoff=(0.5,0.8)):
    ## implementation of clusterJson using csr sparce matrices

    ### cutoff value is in distance (i.e. if two clusters are further than cutoff value, similarity is 0),
    ##  larger cutoff values are more permissive
    bgcs = set() # contains the indices (as floats!) of all the BGCs in this class
    simDict = {}
    # Doing this so it only has to go through the matrix once
    for row in matrix:
        gc1 = row[0]
        gc2 = row[1]
        distance = row[3]
        bgcs.add(gc1)
        bgcs.add(gc2)
        if distance < 1.0:
            similarity = 1 - distance
        else:
            similarity = 0
        gcSimilarities = simDict.setdefault(gc1, {})
        gcSimilarities[gc2] = similarity
    clanClassificationCutoff, clanDistanceCutoff = clanCutoff
    if clusterClans and verbose:
        print('Clustering Clans Enabled with paramters clanClassificationCutoff: {}, clanDistanceCutoff: {}'.format(clanClassificationCutoff,clanDistanceCutoff))
    # if we want to classify by clans make sure that the clanCutoff is included in the cutoffs to do AP in
    if clusterClans and clanClassificationCutoff not in cutoffs:
        cutoffs.append(clanClassificationCutoff)
    # preserve order
    bgcs = sorted(list(bgcs))
    bgc2simIdx = dict(zip(bgcs, range(len(bgcs))))
    pfam_domain_categories = os.path.join(os.path.dirname(os.path.realpath(__file__)), "Pfam-A.clans.tsv")
    pfam_descrs = generatePfamDescriptionsMatrix(pfam_domain_categories)
    for cutoff in cutoffs:
        simMatrix = lil_matrix((len(bgc2simIdx), len(bgc2simIdx)), dtype=np.float32)
        for bgc1 in bgcs:
            # first make sure it is similar to itself
            simMatrix[bgc2simIdx[bgc1],bgc2simIdx[bgc1]] = 1
            for bgc2 in simDict.get(bgc1,{}).keys():
                # you might get 0 values if there were matrix entries under the cutoff don't need to input these in
                # the sparse matrix
                if simDict[bgc1][bgc2] > 1-cutoff:
                    # Ensure symmetry
                    simMatrix[bgc2simIdx[bgc1], bgc2simIdx[bgc2]] = simDict[bgc1][bgc2]
                    simMatrix[bgc2simIdx[bgc2], bgc2simIdx[bgc1]] = simDict[bgc1][bgc2]
        if verbose:
            print("  Clustering (c=" + str(cutoff) + ")")
        labels = pysapc.SAP(damping=damping, max_iter=500,
                            preference='min').fit_predict(simMatrix)
        if verbose:
            print("   ...done")
        numBGCs = len(bgcs)
        bs_distances = [[float("{:.3f}".format(simMatrix[row, col])) for col in range(row+1)] for row in range(numBGCs)]
        #bs_data = [{"id": clusterNames[int(bgc)]} for bgc in bgcs]
        bs_data = []
        bgcJsonDict = {}
        familiesDict = {}
        for idx, label in enumerate(labels):
            members = familiesDict.setdefault(label, [])
            members.append(idx)
            familiesDict[label] = members
            
        ## guarantee order of families and indices for later assignment
        # Families should be indexed from 0 to # families - 1
        familiesDictReIdxd = {idx:members for idx,members in enumerate(familiesDict.values())}
        familyIdxs = range(len(familiesDictReIdxd))
        
        ### Use the 0.5 distance cutoff to cluster clans by default
        if clusterClans and cutoff == clanClassificationCutoff:
            famSimDict = dict()
            for familyI in familyIdxs:
                membersI = familiesDictReIdxd[familyI]
                for familyJ in familyIdxs:
                    membersJ = familiesDictReIdxd[familyJ]
                    famSimDict.setdefault(familyI, {})
                    ### Need to check all of the distances between members, option to get the min, max or mean distances
                    famSimilarities = [0]
                    for memberI in membersI:
                        similarities = [simMatrix[memberI, memberJ] for memberJ in membersJ]
                        # can change this to min or mean
                        famSimilarities.append(sum(similarities, 0.0) / len(similarities))
                    ## can change this as well
                    famSimDict[familyI][familyJ] = sum(famSimilarities,0.0)/len(famSimilarities)

            famSimMatrix = lil_matrix((len(familyIdxs), len(familyIdxs)), dtype=np.float32)
            for familyI in familyIdxs:
                # first make sure it is similar to itself
                famSimMatrix[familyI, familyI] = 1
                for familyJ in famSimDict.get(familyI, {}).keys():
                    # you might get 0 values if there were matrix entries under the cutoff don't need to input these in
                    # the sparse matrix
                    if famSimDict[familyI].get(familyJ, 0) > 1 - clanDistanceCutoff:
                        # Ensure symmetry
                        famSimMatrix[familyI, familyJ] = famSimDict[familyI][familyJ]
                        famSimMatrix[familyJ, familyI] = famSimDict[familyI][familyJ]
            clanLabels = pysapc.SAP(damping=damping, max_iter=500,
                                preference='min').fit_predict(famSimMatrix)
        else:
            clanLabels = []
            
        if len(clanLabels) > 0:
            clansDict = {}
            for idx,clanLabel in enumerate(clanLabels):
                clanMembers = clansDict.setdefault(clanLabel,[])
                clanMembers.append(idx)
                clansDict[clanLabel] = clanMembers
            fam2clan = dict(zip(familyIdxs,clanLabels))
        
        for bgc in bgcs:
            bgcName = clusterNames[int(bgc)]
            bgcJsonDict[bgcName] = {}
            bgcJsonDict[bgcName]["id"] = clusterNames[int(bgc)]
            bgcJsonDict[bgcName]["desc"] = bgc_info[bgcName].description
            bgcJsonDict[bgcName]["start"] = 1
            bgcJsonDict[bgcName]["end"] = bgc_info[bgcName].max_width
            pfdFile = os.path.join(pfd_folder, bgcName + ".pfd")
            fastaFile = os.path.join(bgc_fasta_folder, bgcName + ".fasta")
            orfDict = defaultdict(dict)
            ## read fasta file first to get orfs
            for line in open(fastaFile):
                if line[0] == ">":
                    header = line.strip()[1:].split(':')
                    if header[2]:
                        orfDict[header[0]]["id"] = header[2]
                    elif header[4]:
                        orfDict[header[0]]["id"] = header[4]
                    else:
                        orfDict[header[0]]["id"] = header[0]
                    orfDict[header[0]]["start"] = int(header[6])
                    orfDict[header[0]]["end"] = int(header[7])
                    if header[-1] == '+':
                        orfDict[header[0]]["strand"] = 1
                    else:
                        orfDict[header[0]]["strand"] = -1
                    orfDict[header[0]]["domains"] = []
            ## now read pfd file to add the domains to each of the orfs
            for line in open(pfdFile):
                entry = line.split('\t')
                orf = entry[-1].strip().split(':')[0]
                pfamID = entry[5].split('.')[0]
                pfamDescr = pfam_descrs.get(pfamID,None)
                if pfamDescr:
                    orfDict[orf]["domains"].append({'code': '{} : {}'.format(pfamID,pfamDescr),'start':int(entry[3]),'end':int(entry[4]),'bitscore': float(entry[1])})
                else:
                    orfDict[orf]["domains"].append({'code': entry[5], 'start': int(entry[3]), 'end': int(entry[4]), 'bitscore': float(entry[1])})
            bgcJsonDict[bgcName]['orfs'] = list(orfDict.values())
        bs_data = [bgcJsonDict[clusterNames[int(bgc)]] for bgc in bgcs]
        
        if len(clanLabels) > 0:
            bs_families = [{"id": "FAM_{:03d}".format(family), 'members': members, "clan": "CLAN_{:03d}".format(fam2clan[family])}
                           for family, members in familiesDictReIdxd.items()]
            bs_clans = [{"id": "CLAN_{:03d}".format(clan), 'members': members}
                           for clan, members in enumerate(clansDict.values())]
        else:
            bs_families = [{"id": "FAM_{:03d}".format(family), 'members': members, }
                           for family, members in familiesDictReIdxd.items()]
        
        # column1: BGC, column2: clustering pseudo family
        if verbose:
            print("  Writing clustering file")
        with open("{}_clustering_c{:4.2f}.tsv".format(outputFileBase, cutoff), "w") as clustering_file:
            clustering_file.write('##BGC Name\tFamily Number\n')
            for familyNum,family in enumerate(bs_families):
                for label in family['members']:
                    clustering_file.write('{}\t{}\n'.format(clusterNames[int(bgcs[label])],familyNum))
            #
            # i = 0
            # for label in familiesDictReIdxd:
            #     # moved down so counting starts at 0 (for consistency with everything else)
            #     for x in familiesDictReIdxd[label]:
            #         clustering_file.write(clusterNames[int(bgcs[x])] + "\t" + str(i) + "\n")
            #     i += 1
        if verbose:
            print("  Writing JS file")
        outputFile = "{}_cutoff{:4.2f}.js".format(outputFileBase,cutoff)
        with open(outputFile, 'w') as outfile:
            outfile.write("var bs_similarity={}\n".format(str(bs_distances)))
            outfile.write("var bs_data={}\n".format(str(bs_data)))
            outfile.write("var bs_families={}".format(str(bs_families)))
            if len(clanLabels) > 0:
                outfile.write("\nvar bs_clans={}\n".format(str(bs_clans)))
        if len(clanLabels) > 0:
            if verbose:
                print("   Writing Clans file")
            outputFile = "{}_clans_{:4.2f}_{:4.2f}.tsv".format(outputFileBase,clanClassificationCutoff,clanDistanceCutoff)
            with open(outputFile,'w') as clansFile:
                clansFile.write('##BGC Name\tClan Number\tFamily Number\n')
                for clanNum,clan in enumerate(bs_clans):
                    for familyNum in clan['members']:
                        for BGC_label in bs_families[familyNum]['members']:
                            clansFile.write('{}\t{}\t{}\n'.format(clusterNames[int(bgcs[BGC_label])],clanNum,familyNum))


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
    parser = ArgumentParser()
    
    parser.add_argument("-o", "--outputdir", dest="outputdir", default="", required=True,
                      help="Output directory, this will contain your pfd, pfs, network and hmmscan output files.")
    parser.add_argument("-i", "--inputdir", dest="inputdir", default=os.path.dirname(os.path.realpath(__file__)),
                      help="Input directory of gbk files, if left empty, all gbk files in current and lower directories will be used.")
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(),
                      help="Set the number of cores the script may use (default: use all available cores)")
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="Prints more detailed information. Toggle to true.")
    parser.add_argument("--include_singletons", dest="include_singletons", action="store_true", default=False,
                      help="Include nodes that have no edges to other nodes from the network. Toggle to activate.")
    parser.add_argument("-d", "--domain_overlap_cutoff", dest="domain_overlap_cutoff", default=0.1,
                      help="Specify at which overlap percentage domains are considered to overlap.")
    parser.add_argument("-m", "--min_bgc_size", dest="min_bgc_size", default=0,
                      help="Provide the minimum size of a BGC to be included in the analysis. Default is 0 base pairs")
    
    parser.add_argument("-s", "--samples", dest="samples", action="store_true", default=False, help="Separate the input files into samples according to their containing folder within the input folder. Toggle to activate")
    
    parser.add_argument("--no_all", dest="no_all", action="store_true", default=False, help="By default, BiG-SCAPE uses a single data set comprised of all input files available recursively within the input folder. Toggle to disactivate this behaviour (in that case, if the --samples parameter is not activated, BiG-SCAPE will not create any network file)")
    
    parser.add_argument("--mix", dest="mix", action="store_true", default=False, help="By default, BiG-SCAPE separates the analysis according to the BGC product (PKS Type I, NRPS, RiPPs, etc.) and will create network directories for each class. Toggle to include an analysis mixing all classes")

    parser.add_argument("--cluster_family", dest="cluster_family",action="store_true", default=False, help="BiG-SCAPE will perform a second layer of clustering and attempt to group families assigned from clustering with cutoff of 0.5 to clans")

    parser.add_argument("--clan_cutoff",dest="clan_cutoff",default=[0.5,0.8], type=float, choices=[FloatRange(0.0, 1.0)],nargs=2,
                        help="Cutoff Parameters for which clustering families into clans will be performed in raw distance.\
                             First value is the cutoff value family assignments for BGCs used in clan clustering (default: 0.5). \
                             Second value is the cutoff value for clustering families into clans (default: 0.8). \
                             Average linkage for BGCs in a family is used for distances between families. Example: --clan_cutoff 0.5 0.8)")

    parser.add_argument("--hybrids", dest="hybrids", action="store_true", 
                        default=False, help="Toggle to also add BGCs with hybrid\
                        predicted products from the PKS/NRPS Hybrids and Others\
                        classes to each subclass (e.g. a 'terpene-nrps' BGC from\
                        Others would be added to the Terpene and NRPS classes")
    
    parser.add_argument("--local", dest="local", action="store_true", default=False, help="Activate local mode. BiG-SCAPE will change the logic in the distance calculation phase to try to perform local alignments of shorter, 'fragmented' BGCs by finding the maximum overlap in domain content.")

    parser.add_argument("--no_classify", dest="no_classify", action="store_true", default=False, help="By default, BiG-SCAPE classifies the output files analysis based on the BGC product. Toggle to deactivate (note that if the --mix parameter is not activated, BiG-SCAPE will not create any network file).")

    parser.add_argument("--banned_classes", nargs='+', dest="banned_classes", default=[], choices=["PKSI", "PKSother", "NRPS", "RiPPs", "Saccharides", "Terpene", "PKS-NRP_Hybrids", "Others"], help="Classes that should NOT be included in the classification. E.g. \"--banned_classes PKSI PKSOther\"")

    parser.add_argument("--pfam_dir", dest="pfam_dir",
                      default=os.path.dirname(os.path.realpath(__file__)), 
                      help="Location of hmmpress-processed Pfam files. Default is same location of BiG-SCAPE")
    parser.add_argument("--anchorfile", dest="anchorfile", default=os.path.join(os.path.dirname(os.path.realpath(__file__)),"anchor_domains.txt"),
                      help="Provide a custom location for the anchor domains file, default is anchor_domains.txt.")
    parser.add_argument("--exclude_gbk_str", dest="exclude_gbk_str", default="",
                      help="If this string occurs in the gbk filename, this file will not be used for the analysis.")
    
    parser.add_argument("--mafft_pars", dest="mafft_pars", default="",
                      help="Add single/multiple parameters for MAFFT specific enclosed by quotation marks e.g. \"--nofft --parttree\"")
    parser.add_argument("--al_method", dest="al_method", default="--retree 2",
                      help="alignment method for MAFFT, if there's a space in the method's name, enclose by quotation marks. default: \"--retree 2\" corresponds to the FFT-NS-2 method")
    parser.add_argument("--maxiterate", dest="maxit", default=1000,
                      help="Maxiterate parameter in MAFFT, default is 1000, corresponds to the FFT-NS-2 method")
    parser.add_argument("--mafft_threads", dest="mafft_threads", default=0,
                      help="Set the number of threads in MAFFT, -1 sets the number of threads as the number of physical cores. Default: same as --cores parameter")
    parser.add_argument("--use_mafft", dest="use_mafft", action="store_true", default=False, help="Use MAFFT instead of hmmalign for multiple alignment of domain sequences")
    parser.add_argument("--force_hmmscan", dest="force_hmmscan", action="store_true", default=False, 
                      help="Force domain prediction using hmmscan even if BiG-SCAPE finds processed domtable files (e.g. to use a new version of PFAM).")
    parser.add_argument("--skip_ma", dest="skip_ma", action="store_true", default=False, 
                      help="Skip multiple alignment of domains' sequences. Use if alignments have been generated in a previous run.")
    parser.add_argument("--skip_all", dest="skip_all", action="store_true",
                      default = False, help = "Only generate new network files. ")
    parser.add_argument("--cutoffs", dest="cutoffs", nargs="+", default=[0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80], type=float, choices=[FloatRange(0.0, 1.0)],
                      help="Generate networks using multiple raw distance cutoff values, example: \"0.1, 0.25, 0.5, 1.0\". Default: all values from 0.10 to 0.80 with 0.05 intervals.")

    options = parser.parse_args()
    return options


if __name__=="__main__":
    options = CMD_parser()
    
    class bgc_data:
        def __init__(self, accession_id, description, product, records, max_width, biosynthetic_genes, contig_edge):
            # These two properties come from the genbank file:
            self.accession_id = accession_id
            self.description = description
            # AntiSMASH predicted class of compound:
            self.product = product
            # number of records in the genbank file (think of multi-locus BGCs):
            self.records = records
            # length of largest record (it will be used for ArrowerSVG):
            self.max_width = int(max_width)
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
    global verbose
    global BGCs
    
    # contains the type of the final product of the BGC (as predicted by AntiSMASH), 
    # as well as the definition line from the BGC file. Used in the final network files.
    global output_folder

    global pfam_dir
    global timings_file
    global cores
    global local
    global clusterNames, bgcClassNames
    
    include_singletons = options.include_singletons
    
    cores = int(options.cores)
    
    options_all = not options.no_all
    options_samples = options.samples
    
    options_mix = options.mix
    options_classify = not options.no_classify
    local = options.local
    
    cutoff_list = options.cutoffs
    for c in cutoff_list:
        if c <= 0.0 or c > 1.0:
            cutoff_list.remove(c)

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
    
    if options.mafft_threads == 0:
        options.mafft_threads = options.cores
    else:
        options.mafft_threads = int(options.mafft_threads)
                    
    verbose = options.verbose
    
    networks_folder_all = "networks_all"
    networks_folder_samples = "networks_samples"    
    if options.hybrids:
        networks_folder_all += "_hybrids"
        networks_folder_samples += "_hybrids"
    if local:
        networks_folder_all += "_local"
        networks_folder_samples += "_local"
    
    if options.skip_all and options.skip_ma:
        print("Overriding --skip_ma with --skip_all parameter")
        options.skip_hmmscan = False
        options.skip_ma = False
    
    time1 = time.time()
    
    # Make the following available for possibly deleting entries within parseHmmScan
    global genbankDict, gbk_files, sampleDict, clusters, baseNames
    
    
    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    print("\n\n   - - Processing input files - -")
    
    create_directory(output_folder, "Output", False)
    write_parameters(output_folder, options)
    bgc_fasta_folder = os.path.join(output_folder, "fasta")
    create_directory(bgc_fasta_folder, "BGC fastas", False)
    
    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    bgc_info = {} # Stores, per BGC: predicted type, gbk Description, number of records, width of longest record, GenBank's accession
    genbankDict = get_gbk_files(options.inputdir, bgc_fasta_folder, int(options.min_bgc_size), options.exclude_gbk_str, bgc_info)

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
    
    domtable_folder = os.path.join(output_folder, "domtable")
    pfs_folder = os.path.join(output_folder, "pfs")
    pfd_folder = os.path.join(output_folder, "pfd")    
    domains_folder = os.path.join(output_folder, "domains")
    svg_folder = os.path.join(output_folder, "SVG")
    
    create_directory(domtable_folder, "Domtable", False)
    create_directory(domains_folder, "Domains", False)
    create_directory(pfs_folder, "pfs", False)
    create_directory(pfd_folder, "pfd", False)
    create_directory(svg_folder, "SVG", False)

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
            with open(os.path.join(output_folder, "BGCs.dict"), "r") as BGC_file:
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
        with open(os.path.join(output_folder, "BGCs.dict"), "wb") as BGC_file:
            pickle.dump(BGCs, BGC_file)
            BGC_file.close()

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
        
        print("  Parsing hmm file for domain names")
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
        fasta_domains = set(glob(os.path.join(domains_folder,"*.fasta")))
        
        # compare with .algn set of files. Maybe resuming is possible if
        # no new sequences were added
        if try_MA_resume:
            temp_aligned = set(glob(os.path.join(domains_folder, "*.algn")))
            
            if len(temp_aligned) > 0:
                print(" Found domain fasta files without corresponding alignments")
                
                for a in temp_aligned:
                    if os.path.getsize(a) > 0:
                        fasta_domains.remove(a[:-5]+".fasta")
            
            temp_aligned.clear()
        
        # Try to further reduce the set of domain fastas that need alignment
        sequence_tag_list = set()
        header_list = []
        fasta_domains_temp = fasta_domains.copy()
        for domain_file in fasta_domains_temp:
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
                fasta_domains.remove(domain_file)
                if verbose:
                    print(" Skipping Multiple Alignment for {} (appears only in one BGC)".format(domain_name))
        
        sequence_tag_list.clear()
        del header_list[:]
        
        fasta_domains_temp.clear()
            
        # Do the multiple alignment
        if len(fasta_domains) > 0:
            if not options.use_mafft:
                print("\n Using hmmalign")
                launch_hmmalign(cores, fasta_domains)
                                        
            else:
                print("\n Using MAFFT")
                for domain in fasta_domains:
                    domain_name = domain[:-6]
                    # Multiple alignment of all domain sequences
                    run_mafft(options.al_method, options.maxit, options.mafft_threads, options.mafft_pars, domain_name)
    
            # verify all tasks were completed by checking existance of alignment files
            for domain in fasta_domains:
                if not os.path.isfile(domain[:-6]+".algn"):
                    print("   WARNING, {}.algn could not be found (possible issue with aligner).".format(domain[:-6]))
                       
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

    clusterNames = tuple(sorted(list(clusters)))
    clusterNames2idx = dict(zip(clusterNames,range(len(clusterNames))))

    # Try to make default analysis using all files found inside the input folder
    if options_all:
        print("\nGenerating distance network files with ALL available input files")
    
        # create output directory
        create_directory(os.path.join(output_folder, networks_folder_all), "Networks_all", False)
    
        # Making network files mixing all classes
        if options_mix:
            print("\n Mixing all BGC classes")
            
            # only choose from valid classes
            mix_set = []
            for clusterIdx,cluster in enumerate(clusterNames):
                product = bgc_info[cluster].product
                predicted_class = sort_bgc(product)
                if predicted_class.lower() in valid_classes:
                    mix_set.append(clusterIdx)
            
            # Create an additional file with the list of all clusters in the class + other info
            print("   Writing annotation file")
            path_list = os.path.join(output_folder, networks_folder_all, "Network_Annotations_ALL_mix.tsv")
            with open(path_list, "w") as list_file:
                list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\n")
                for idx in mix_set:
                    bgc = clusterNames[idx]
                    product = bgc_info[bgc].product
                    list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product)]) + "\n")
            
            print("  Calculating all pairwise distances")
            pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
            del mix_set[:]
            cluster_pairs = [(x, y, -1) for (x, y) in pairs]
            pairs.clear()
            network_matrix_mix = generate_network(cluster_pairs, cores)
            del cluster_pairs[:]
                
            print("  Writing output files")
            pathBase = os.path.join(output_folder, networks_folder_all, "all_mix")
            filenames = []
            for cutoff in cutoff_list:
                filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
            clusterJsonBatch(pathBase, network_matrix_mix, cutoffs=cutoff_list,clusterClans=options.cluster_family,clanCutoff=options.clan_cutoff)
            cutoffs_and_filenames = list(zip(cutoff_list, filenames))
            del filenames[:]
            write_network_matrix(network_matrix_mix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
            
            del network_matrix_mix[:]
            
        # Making network files separating by BGC class
        if options_classify:
            print("\n Working for each BGC class")
            
            # reinitialize BGC_classes to make sure the bgc lists are empty
            BGC_classes = defaultdict(list)
        
            # Preparing gene cluster classes
            print("  Sorting the input BGCs\n")
            for clusterIdx,cluster in enumerate(clusterNames):
                product = bgc_info[cluster].product
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
                            
                        for subclass in subclasses:
                            BGC_classes[subclass].append(clusterIdx)
                        subclasses.clear()

            # only make folders for the BGC_classes that are found
            for bgc_class in BGC_classes:
                folder_name = bgc_class
                    
                print("\n  {} ({} BGCs)".format(folder_name, str(len(BGC_classes[bgc_class]))))
                
                # create output directory   
                create_directory(os.path.join(output_folder, networks_folder_all, folder_name), "  All - " + bgc_class, False)
                
                # Create an additional file with the final list of all clusters in the class
                print("   Writing annotation files")
                path_list = os.path.join(output_folder, networks_folder_all, folder_name, "Network_Annotations_All_" + folder_name + ".tsv")
                with open(path_list, "w") as list_file:
                    list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\n")
                    for idx in BGC_classes[bgc_class]:
                        bgc = clusterNames[idx]
                        product = bgc_info[bgc].product
                        list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product)]) + "\n")
            
                    
                if len(BGC_classes[bgc_class]) > 1:
                    print("   Calculating all pairwise distances")
                    pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                    del BGC_classes[bgc_class][:]
                    cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                    pairs.clear()
                    network_matrix = generate_network(cluster_pairs, cores)
                    del cluster_pairs[:]
                        
                    print("   Writing output files")
                    pathBase = os.path.join(output_folder, networks_folder_all, folder_name, "all" + folder_name)
                    filenames = []
                    for cutoff in cutoff_list:
                        filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
                    cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                    del filenames[:]
                    clusterJsonBatch(pathBase, network_matrix, cutoffs=cutoff_list,clusterClans=options.cluster_family,clanCutoff=options.clan_cutoff)
                    write_network_matrix(network_matrix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
                    
                    del network_matrix[:]
                

    # Try to make analysis for each sample
    if options_samples:
        network_matrix_sample = []
        
        if len(sampleDict) == 1 and options_all:
            print("\nNOT generating networks per sample (only one sample, covered in the all-vs-all case)")
        else:
            print("\nGenerating distance network files for each sample")
            
            # create output directory for all samples
            create_directory(os.path.join(output_folder, networks_folder_samples), "Samples", False)
            
            for sample, sampleClusters in sampleDict.items():
                print("\n Sample: " + sample)
                if len(sampleClusters) == 1:
                    print(" Warning: Sample size = 1 detected. Not generating network for this sample ({})".format(sample))
                else:
                    # create output directory for this sample
                    create_directory(os.path.join(output_folder, networks_folder_samples, sample), " Samples - " + sample, False)
                        
                    # Making network files mixing all classes
                    if options_mix:
                        print("\n  Mixing all BGC classes")
                        
                        mix_set = []
                        for clusterIdx, sampleCluster in enumerate(sampleClusters):
                            product = bgc_info[sampleCluster].product
                            predicted_class = sort_bgc(product)
                            if predicted_class.lower() in valid_classes:
                                mix_set.append(clusterIdx)
                        
                        # Create an additional file with the list of all clusters in the class + other info
                        print("   Writing annotation files")
                        path_list = os.path.join(output_folder, networks_folder_samples, sample, "Network_Annotations_Sample_" + sample + "_mix.tsv")
                        with open(path_list, "w") as list_file:
                            list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\n")
                            for idx in mix_set:
                                bgc = clusterNames[idx]
                                product = bgc_info[bgc].product
                                list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product)]) + "\n")
            
                        pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
                        del mix_set[:]
                        cluster_pairs = [(x, y, -1) for (x, y) in pairs]
                        pairs.clear()
                        network_matrix_sample = generate_network(cluster_pairs, cores)
                        del cluster_pairs[:]

                        print("   Writing output files")
                        pathBase = os.path.join(output_folder, networks_folder_samples, sample, "sample_" + sample + "_mix")
                        filenames = []
                        for cutoff in cutoff_list:
                            filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
                        cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                        clusterJsonBatch(pathBase, network_matrix_sample, cutoffs=cutoff_list,clusterClans=options.cluster_family,clanCutoff=options.clan_cutoff)
                        write_network_matrix(network_matrix_sample, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
                        
                        del network_matrix_sample[:]
                    
                    # Making network files separating by BGC class
                    if options_classify:
                        print("\n  Working for each BGC class")

                        # reinitialize BGC_classes to make sure the bgc lists are empty
                        BGC_classes = defaultdict(list)
                    
                        # Preparing gene cluster classes
                        print("   Sorting the input BGCs\n")
                        for cluster in sampleClusters:
                            product = bgc_info[cluster].product
                            predicted_class = sort_bgc(product)
                            if predicted_class.lower() in valid_classes:
                                BGC_classes[predicted_class].append(clusterNames2idx[cluster])
                            
                            # possibly add hybrids to 'pure' classes
                            if options.hybrids:
                                if predicted_class == "PKS-NRP_Hybrids":
                                    if "nrps" in valid_classes:
                                        BGC_classes["NRPS"].append(clusterNames2idx[cluster])
                                    if "t1pks" in product and "pksi" in valid_classes:
                                        BGC_classes["PKSI"].append(clusterNames2idx[cluster])
                                    if "t1pks" not in product and "pksother" in valid_classes:
                                        BGC_classes["PKSother"].append(clusterNames2idx[cluster])
                                
                                if predicted_class == "Others" and "-" in product:
                                    subclasses = set()
                                    for subproduct in product.split("-"):
                                        subclass = sort_bgc(subproduct)
                                        if subclass.lower() in valid_classes:
                                            subclasses.add(subclass)
                                        
                                    for subclass in subclasses:
                                        BGC_classes[subclass].append(clusterNames2idx[cluster])
                                    subclasses.clear()

                        for bgc_class in BGC_classes:
                            folder_name = bgc_class
                                
                            print("\n   " + folder_name + " (" + str(len(BGC_classes[bgc_class])) + " BGCs)")
                            network_matrix_sample = []
                            
                            # create output directory
                            create_directory(os.path.join(output_folder, networks_folder_samples, sample, folder_name), "   Sample " + sample + " - " + bgc_class, False)

                            # Create an additional file with the final list of all clusters in the class
                            path_list = os.path.join(output_folder, networks_folder_samples, sample, folder_name, "Network_Annotations_Sample_" + sample + "_" + folder_name + ".tsv")
                            with open(path_list, "w") as list_file:
                                list_file.write("BGC\tAccesion ID\tDescription\tProduct Prediction\tBiG-SCAPE class\n")
                                for idx in BGC_classes[bgc_class]:
                                    bgc = clusterNames[idx]
                                    product = bgc_info[bgc].product
                                    list_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product)]) + "\n")

                            if len(BGC_classes[bgc_class]) > 1:
                                pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                                del BGC_classes[bgc_class][:]
                                cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                                pairs.clear()
                                network_matrix_sample = generate_network(cluster_pairs, cores)
                                del cluster_pairs[:]
                                print("    Writing output files")
                                pathBase = os.path.join(output_folder, networks_folder_samples, sample, folder_name,
                                                        "sample_" + sample + "_" + folder_name)
                                filenames = []
                                for cutoff in cutoff_list:
                                    filenames.append("{}_c{:.2f}.network".format(pathBase, cutoff))
                                cutoffs_and_filenames = list(zip(cutoff_list, filenames))
                                clusterJsonBatch(pathBase, network_matrix_sample, cutoffs=cutoff_list,clusterClans=options.cluster_family,clanCutoff=options.clan_cutoff)
                                write_network_matrix(network_matrix_sample, cutoffs_and_filenames, include_singletons, clusterNames,bgc_info)
                                del network_matrix_sample[:]
                                
                            del BGC_classes[bgc_class][:]

    pickle.dump(bgc_info,open(os.path.join(output_folder,'bgc_info.dict'),'wb'))
    runtime = time.time()-time1
    runtime_string = "\n\n\tMain function took {:.3f} s".format(runtime)
    with open(os.path.join(output_folder, "runtimes.txt"), 'a') as timings_file:
        timings_file.write(runtime_string + "\n")
    print(runtime_string)
    
