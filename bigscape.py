#!/usr/bin/env python


"""
BiG-SCAPE

PI: Marnix Medema

Developers:
Marley Yeong                    marleyyeong@live.nl
Jorge Navarro
Emmanuel (Emzo) de los Santos

Dependencies: hmmer, biopython, mafft, munkres

Usage:   Please see `python bigscape.py -h`

Example: python bigscape.py -c 8 --pfam_dir ./ -i ./inputfiles -o ./results

Status: development/testing

"""

import cPickle as pickle  # for storing and retrieving dictionaries
from math import exp, log
import os
import subprocess
import sys
import time
from glob import glob
from itertools import combinations
from multiprocessing import Pool, cpu_count
from optparse import OptionParser

from Bio import SeqIO
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix

from functions import *
from munkres import Munkres

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


def get_gbk_files(inputdir, min_bgc_size, exclude_gbk_str, gbk_group):
    """Searches given directory for genbank files recursively, will assume that
    the genbank files that have the same name are the same genbank file. 
    Returns a dictionary that contains the names of the clusters found as keys
    and a list that contains [0] a path to the genbank file and [1] the 
    samples that the genbank file is a part of.
    return: {cluster_name:[genbank_path,[s_a,s_b...]]}
    """

    genbankDict = {}

    file_counter = 0

    print("\nImporting GenBank files")
    if exclude_gbk_str != "":
        print(" Skipping files with '" + exclude_gbk_str + "' in their filename")

    # this doesn't seem to make a difference. Confirm
    # if inputdir != "" and inputdir[-1] != "/":
    # inputdir += "/"
    current_dir = ""
    for dirpath, dirnames, filenames in os.walk(inputdir):
        head, tail = os.path.split(dirpath)
        if current_dir != tail:
            current_dir = tail

        genbankfilelist = []

        # avoid double slashes
        dirpath_ = dirpath[:-1] if dirpath[-1] == os.sep else dirpath

        for fname in filenames:
            fname_parse = fname.split('.')
            
            if fname_parse[-1] != "gbk":
                continue
            
            clusterName = '.'.join(fname_parse[:-1])
            
            if exclude_gbk_str != "" and exclude_gbk_str in fname:
                print(" Skipping file " + fname)
                continue
            
            if " " in fname:
                sys.exit("\nError: Input GenBank files should not have spaces in their filenames as HMMscan cannot process them properly ('too many arguments').")
                
            
            with open(os.path.join(dirpath_, fname), "r") as f:
                try:
                    # basic file verification. Substitutes check_data_integrity
                    # look into SeqIO.parse for multiple records [FUTURE]
                    record = SeqIO.read(f, "genbank")
                except ValueError as e:
                    print("   Error with file " + os.path.join(dirpath_, fname) + ": \n    '" + str(e) + "'")
                    print("    (This file will be excluded from the analysis)")
                    continue
                else:
                    group = "no type"
                    for feature in record.features:
                        if "cluster" in feature.type and "product" in feature.qualifiers:
                            group = ",".join(feature.qualifiers["product"])
                    gbk_group[clusterName] = (group, record.description)
                    
                    bgc_size = len(record.seq)
                    if bgc_size > min_bgc_size:  # exclude the bgc if it's too small
                        file_counter += 1
                        
                        if clusterName in genbankDict.keys():
                            # current_dir gets to be the name of the sample
                            genbankDict[clusterName][1].add(current_dir) 
                        else:
                            # location of first instance of the file is genbankDict[clustername][0]
                            genbankDict.setdefault(clusterName, [os.path.join(dirpath_, fname), set([current_dir])])
                            
                        if verbose == True:
                            print("  Adding " + fname + " (" + str(bgc_size) + " bps)")
                    else:
                        print(" Discarding " + clusterName +  " (size less than " + str(min_bgc_size) + " bp, was " + str(bgc_size) + ")")
                        
            # else: The file does not have the gbk extension. Skip it
    
    if file_counter == 0:
        sys.exit("\nError: There are no files to process")
        
    if file_counter == 1:
        sys.exit("\nError: Only one file found. Please input at least two files")
    
    if verbose:
        print("\n Starting with " + str(file_counter) + " files")

    return genbankDict

def timeit(f):
    def wrap(*args):
        insignificant_runtime = 1 #prevents an overload 
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        runtime = time2-time1
        
        runtime_string = '\t%s function took %0.3f s' % (f.func_name, runtime)
        
        if runtime > insignificant_runtime:
            timings_file.write(runtime_string + "\n")
            print runtime_string
            
        return ret
    return wrap


def generate_dist_matrix(parms):    
    #Get the values from the parameters
    cluster1 = parms[0]
    cluster2 = parms[1]
    dist_method = parms[2]
    anchor_domains = parms[3]
    
    try:
        domain_list_A = DomainList[cluster1]
        domain_list_B = DomainList[cluster2]
    except KeyError:
        if verbose:
            print(" Warning: domain list for " + cluster1 + " or " + cluster2 + " was not found. Extracting from pfs files")
        
        cluster_file1 = os.path.join(output_folder, cluster1 + ".pfs")
        cluster_file2 = os.path.join(output_folder, cluster2 + ".pfs")
        
        domain_list_A = get_domain_list(cluster_file1)
        domain_list_B = get_domain_list(cluster_file2)
    
    # this really shouldn't happen if we've filtered domain-less gene clusters already
    if len(domain_list_A) == 0 or len(domain_list_B) == 0:
        print("   Warning: Regarding distance between clusters " + cluster1 + " and " + cluster2 + ":")
        if len(domain_list_A) == 0 and len(domain_list_B) == 0:
            print("   None have identified domains. Distance cannot be calculated")
        elif (domain_list_A) == 0:            
            print("   Cluster " + cluster1 + " has no identified domains. Distance set to 1")
        else:
            print("   Cluster " + cluster2 + " has no identified domains. Distance set to 1")
            
        # last two values (S, Sa) should really be zero but this could give rise to errors when parsing 
        # the network file (unless we catched the case S = Sa = 0
        return [cluster1, cluster2, group_dct[cluster1][0], group_dct[cluster1][1], \
            group_dct[cluster2][0],group_dct[cluster2][1], '0.0', '1.0', '0.0', '0.0', '0.0', '0.0', "1.0", "1.0", "1", "1"] 
    

    if dist_method == "seqdist":
        dist, jaccard, dds, gk, rDDSna, rDDS, S, Sa = cluster_distance(cluster1, cluster2, domain_list_A, domain_list_B, anchor_domains) #sequence dist
    elif dist_method == "domain_dist":
        dist, jaccard, dds, gk, rDDSna, rDDS, S, Sa = Distance_modified(domain_list_A, domain_list_B, 0, 4) #domain dist
        
    if dist == 0:
        logscore = float("inf")
    else:
        logscore = 0
        try:
            logscore = -log(dist, 2) #Write exception, ValueError
        except ValueError:
            print "calculating the logscore with distance", dist, "failed in function generate_network, using distance method", dist_method, networkfilename
            
    #clustername1 clustername2 group1, group2, -log2score, dist, squared similarity, j, dds, gk
    network_row = [str(cluster1), str(cluster2), group_dct[cluster1][0],group_dct[cluster1][1], \
        group_dct[cluster2][0],group_dct[cluster2][1], str(logscore), str(dist), str((1-dist)**2), \
        jaccard, dds, gk, rDDSna, rDDS, S, Sa]
    
    return network_row
    
        
@timeit
def generate_network(cluster_pairs, cores):
    #Contents of the network file: clustername1 clustername2, group1, group2, -log2score, dist, squared similarity
    "saves the distances as the log2 of the similarity"
    pool = Pool(cores, maxtasksperchild=500)
    
    #Assigns the data to the different workers and pools the results back into
    # the network_matrix variable
    network_matrix = pool.map(generate_dist_matrix, cluster_pairs)
    
    # --- Serialized version of distance calculation ---
    # For the time being, use this if you have memory issues
    #network_matrix = []
    #for pair in cluster_pairs:
      #network_matrix.append(generate_dist_matrix(pair))

    # use a dictionary to store results
    network_matrix_dict = {}
    for row in network_matrix:
        network_matrix_dict[row[0], row[1]] = row[2:]
    
    return network_matrix_dict
    
    
       
def cluster_distance(A, B, A_domlist, B_domlist, anchor_domains): 
    """Compare two clusters using information on their domains, and the sequences of the domains"""    

    temp_domain_fastas = {}

    intersect = set(A_domlist).intersection(B_domlist)
    not_intersect = set(A_domlist).symmetric_difference(set(B_domlist))
    
    
    # JACCARD INDEX
    Jaccard = len(intersect)/ float( len(set(A_domlist)) + len(set(B_domlist)) - len(intersect))


    # DDS INDEX
    #domain_difference: Difference in sequence per domain. If one cluster doesn't have a domain at all, but the other does, 
    #this is a sequence difference of 1. If both clusters contain the domain once, and the sequence is the same, there is a seq diff of 0.
    #S: Max occurence of each domain
    domain_difference_anchor,S_anchor = 0,0 
    domain_difference,S = 0,0 
    
    # Case 1
    for unshared_domain in not_intersect: #no need to look at seq identity, since these domains are unshared
        #for each occurence of an unshared domain do domain_difference += count of domain and S += count of domain
        unshared_occurrences = []
        try:
            unshared_occurrences = BGCs[A][unshared_domain]
        except KeyError:
            unshared_occurrences = BGCs[B][unshared_domain]
            
        # don't look at domain version, hence the split
        if unshared_domain.split(".")[0] in anchor_domains:
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
        
        num_copies_a = len(specific_domain_list_A)
        num_copies_b = len(specific_domain_list_B)
        
        temp_domain_fastas.clear()
        
        # Case 2: The shared domains occurs only once in each gene cluster
        #if len(specific_domain_list_A+specific_domain_list_B) == 2: #The domain occurs only once in both clusters
        #   print(this is case 2)
        # Case 3: The domain occurs more than once in both clusters
        #else:
        #   print(this is case 3)
        
        accumulated_distance = 0
            
        # Fill distance matrix between domain's A and B versions
        DistanceMatrix = [[1 for col in range(num_copies_b)] for row in range(num_copies_a)]
        
        for domsa in range(num_copies_a):
            for domsb in range(num_copies_b):
                sequence_tag_a = specific_domain_list_A[domsa]
                sequence_tag_b = specific_domain_list_B[domsb]
                
                seq_length = 0
                matches = 0
                gaps = 0
                
                try:
                    aligned_seqA = AlignedDomainSequences[sequence_tag_a]
                    aligned_seqB = AlignedDomainSequences[sequence_tag_b]
                    
                except KeyError:
                    # For some reason we don't have the multiple alignment from MAFFT. 
                    # Try manual alignment
                    if shared_domain not in missing_aligned_domain_files and verbose:
                        # this will print everytime an unfound <domain>.algn is not found for every
                        # distance calculation (but at least, not for every domain pair!)
                        print("  Warning: " + shared_domain + ".algn not found. Trying pairwise alignment...")
                        missing_aligned_domain_files.append(shared_domain)
                    
                    try:
                        unaligned_seqA = temp_domain_fastas[sequence_tag_a]
                        unaligned_seqB = temp_domain_fastas[sequence_tag_b]
                    except KeyError:
                        # parse the file for the first time and load all the sequences
                        with open(os.path.join(output_folder, domainsout, shared_domain + ".fasta"),"r") as domain_fasta_handle:
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
                    print("\tWARNING: mismatch in sequences' lengths while calculating sequence identity (" + shared_domain + ")")
                    if verbose:
                        print("\t  Specific domain 1: " + aligned_seqA + " len: " + str(len(seq1)))
                        print("\t  Specific domain 2: " + aligned_seqB + " len: " + str(len(seq2)))
                    seq_length = min(len(aligned_seqA), len(aligned_seqB))
                else:
                    seq_length = len(aligned_seqA)
                    
                for position in range(seq_length):
                    if aligned_seqA[position] == aligned_seqB[position]:
                        if aligned_seqA[position] != "-":
                            matches += 1
                        else:
                            gaps += 1
                            
                DistanceMatrix[domsa][domsb] = 1 - ( float(matches)/float(seq_length-gaps) )
                            
            
        #Only use the best scoring pairs
        Hungarian = Munkres()
        #print "DistanceMatrix", DistanceMatrix
        BestIndexes = Hungarian.compute(DistanceMatrix)
        #print "BestIndexes", BestIndexes
        accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
        #print "accumulated_distance", accumulated_distance
        
        # the difference in number of domains accounts for the "lost" (or not duplicated) domains
        sum_seq_dist = (abs(num_copies_a-num_copies_b) + accumulated_distance)  #essentially 1-sim
        
        if shared_domain.split(".")[0] in anchor_domains: 
            S_anchor += max(num_copies_a,num_copies_b)
            domain_difference_anchor += sum_seq_dist
        else:
            S += max(num_copies_a, num_copies_b)
            domain_difference += sum_seq_dist
        
        
    
    if S_anchor != 0 and S != 0:
        DDS_non_anchor = domain_difference / float(S)
        DDS_anchor = domain_difference_anchor / float(S_anchor)
        
        # Calculate proper, proportional weight to each kind of domain
        non_anchor_prct = S / float(S + S_anchor)
        anchor_prct = S_anchor / float(S + S_anchor)
        
        # boost anchor subcomponent and re-normalize
        non_anchor_weight = non_anchor_prct / (anchor_prct*anchorweight + non_anchor_prct)
        anchor_weight = anchor_prct*anchorweight / (anchor_prct*anchorweight + non_anchor_prct)

        # Use anchorweight parameter to boost percieved rDDS_anchor
        DDS = (non_anchor_weight*DDS_non_anchor) + (anchor_weight*DDS_anchor)
        
    elif S_anchor == 0:
        DDS_non_anchor = domain_difference / float(S)
        DDS_anchor = 0.0
        
        DDS = DDS_non_anchor
        
    else: #only anchor domains were found
        DDS_non_anchor = 0.0
        DDS_anchor = domain_difference_anchor / float(S_anchor)
        
        DDS = DDS_anchor
 
    DDS = 1-DDS #transform into similarity
 
 
    # GK INDEX
    #  calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A_domlist]
    Ar.reverse()
    GK = max([calculate_GK(A_domlist, B_domlist, nbhood), calculate_GK(Ar, B_domlist, nbhood)])
    
    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK) 
    if Distance < 0:
        print("Negative distance detected!")
        print("J: " + str(Jaccard) + "\tDDS: " + str(DDS) + "\tGK: " + str(GK))
        print("Jw: " + str(Jaccardw) + "\tDDSw: " + str(DDSw) + "\tGKw: " + str(GKw))
        sys.exit()
        
    return Distance, Jaccard, DDS, GK, DDS_non_anchor, DDS_anchor, S, S_anchor



@timeit
def run_mafft(al_method, maxit, cores, mafft_pars, domain):
    """Runs mafft program with the provided parameters.
    The specific parameters used to run mafft with are actually very important for the final result.
    Using mafft with the most progressive parameters does indeed affect the quality of the distance.
    It is better to just use the domain information for the distance if more computationally intensive options
    for mafft cannot be used. Setting maxiterate to anything higher than 10 did not make a difference in accuracy in the nrps testset"""
    
    alignment_file = domain + ".algn"
    
    
    mafft_cmd_list = []
    mafft_cmd_list.append("mafft --distout") #distout will save the distance matrix in a ".hat2" file
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


@timeit
def calculate_GK(A, B, nbhood):
    """Goodman and Kruskal's gamma is a measure of rank correlation, i.e., 
    the similarity of the orderings of the data when ranked by each of the quantities."""
    GK = 0.
    if len(set(A) & set(B)) > 1:
        pairsA = set( [(A[i],A[j]) for i in xrange(len(A)-1) for j in xrange(i+1,(i+nbhood if i+nbhood < len(A) else len(A)))] )
        pairsB = set( [(B[i],B[j]) for i in xrange(len(B)-1) for j in xrange(i+1,(i+nbhood if i+nbhood < len(B) else len(B)))] )
        allPairs = set(list(pairsA) + list(pairsB))
        Ns, Nr = 0.,0.
        for p in allPairs:
            if p in pairsA and p in pairsB: Ns += 1
            elif p in pairsA and tuple(p[::-1]) in pairsB: Nr += 1
            elif tuple(p[::-1]) in pairsA and p in pairsB: Nr += 1
            else: pass
        
        if (Nr + Ns) == 0: # this could happen if e.g. only two domains are shared but are farther than nbhood
            gamma = 0
        else:
            gamma = (Ns-Nr) / (Nr+Ns)
        GK = (1+gamma)/2.
    return GK


@timeit
def Distance_modified(clusterA, clusterB, repeat=0, nbhood=4):
    "Modified to work better for 'OBU' detection"
    "Original DDS formula from Lin, Zhu and Zhang (2006)"

    repeats = []

    # delete short and frequent domains
    if repeat==1:
        A = [i for i in clusterA] 
        B = [j for j in clusterB] 
    elif repeat==0:
        A = [i for i in clusterA if i not in repeats]
        B = [j for j in clusterB if j not in repeats]

    if len(A)==0 or len(B)==0: return 1.

    # calculate Jaccard index, modified not to give problems with size differences between clusters
    Jaccard = len(set(A) & set(B)) / float( 2 * min([len(set(A)),len(set(B))]) - len(set(A) & set(B)) )
    #Jaccard = len(set(A) & set(B)) / float( len(set(A)) + len(set(B)) - len(set(A) & set(B)) )

    # calculate domain duplication index
    DDS = 0 #The difference in abundance of the domains per cluster
    S = 0 #Max occurence of each domain
    for p in set(A+B):
        DDS += abs(A.count(p)-B.count(p))
        S += max(A.count(p),B.count(p))
    DDS /= float(S) 
    DDS = exp(-DDS) #transforms the DDS to a value between 0 - 1

    # calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A]
    Ar.reverse()
    GK = max([calculate_GK(A, B, nbhood), calculate_GK(Ar, B, nbhood)]) #100% dissimilarity results in a score of 0.5


    # calculate the distance
    #print "Jaccard", Jaccard
    #print "DDS", DDS
    #print "GK", GK
    Distance = 1 - Jaccardw*Jaccard - DDSw*DDS - GKw*GK
    
    if Distance < 0:
        Distance = 0

    return Distance, Jaccard, DDS, GK

def generateFasta(gbkfilePath,outputdir):
    ## first parse the genbankfile and generate the fasta file for input into hmmscan ##
    outputbase  = gbkfilePath.split(os.sep)[-1].replace(".gbk","")
    if verbose:
        print "   Generating fasta for: ", outputbase
    outputfile = os.path.join(outputdir,outputbase + '.fasta')

    with open(gbkfilePath,"r") as genbankHandle:
        genbankEntry = SeqIO.read(genbankHandle,"genbank")
        CDS_List = (feature for feature in genbankEntry.features if feature.type == 'CDS')

        cds_ctr = 0
        # parse through the CDS lists to make the fasta file for hmmscan, if translation isn't available attempt manual translation
        for CDS in CDS_List:
            cds_ctr += 1
            gene_id = CDS.qualifiers.get('gene',"")
            protein_id = CDS.qualifiers.get('protein_id',"")
            gene_start = max(0,CDS.location.nofuzzy_start)
            gene_end = max(0,CDS.location.nofuzzy_end)
            direction = CDS.location.strand

            if direction == 1:
                strand = '+'
            else:
                strand = '-'

            if 'translation' in CDS.qualifiers.keys():
                prot_seq = CDS.qualifiers['translation'][0]
            # If translation isn't available translate manually, this will take longer
            else:
                genbank_seq = CDS.location.extract(genbank_entry)


                nt_seq = genbank_seq.seq
                if direction == 1:
                    # for protein sequence if it is at the start of the entry assume that end of sequence is in frame
                    # if it is at the end of the genbank entry assume that the start of the sequence is in frame
                    if gene_start == 0:
                        if len(nt_seq) % 3 == 0:
                            prot_seq = nt_seq.translate()
                        elif len(nt_seq) % 3 == 1:
                            prot_seq = nt_seq[1:].translate()
                        else:
                            prot_seq = nt_seq[2:].translate()
                    else:
                        prot_seq = nt_seq.translate()
                # reverse direction
                else:
                    #same logic reverse direction
                    if gene_start == 0:
                        if len(nt_seq) % 3 == 0:
                            prot_seq = nt_seq.translate()
                        elif len(nt_seq) % 3 == 1:
                            prot_seq = nt_seq[:-1].translate()
                        else:
                            prot_seq = nt_seq[:-2].translate()
                    else:
                        prot_seq = nt_seq.translate()

            # write fasta file
            with open(outputfile,'ab') as fastaHandle:
                # final check to see if string is empty
                if prot_seq:
                    fasta_header = outputbase + "_ORF" + str(cds_ctr)+ ":gid:" + str(gene_id) + ":pid:" + str(protein_id) + \
                                   ":loc:" + str(gene_start) + ":" + str(gene_end) + ":strand:" + strand
                    fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this
                    fasta_header = ">"+(fasta_header.replace(" ", "")) #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header
                    fastaHandle.write('%s\n' % fasta_header)
                    fastaHandle.write('%s\n' % prot_seq)
    return outputfile

def runHmmScan(fastaPath, hmmPath, outputdir, verbose):
    ## will run hmmscan command on a fasta file with a single core and generate a domtable file
    hmmFile = os.path.join(hmmPath,"Pfam-A.hmm")
    if os.path.isfile(fastaPath):
        name = fastaPath.split(os.sep)[-1].replace(".fasta","")
        outputName = os.path.join(outputdir, name+".domtable")
        
        hmmscan_cmd = "hmmscan --cpu 1 --domtblout %s --cut_tc %s %s" % (outputName,hmmFile,fastaPath)
        if verbose == True:
            print("  Processing " + name)
            print("   " + hmmscan_cmd)
        subprocess.check_output(hmmscan_cmd, shell=True)

    else:
        sys.exit("Error running hmmscan: Fasta file " + fastaPath + " doesn't exist")

def parseHmmScan(hmmscanResults,outputdir,overlapCutoff):
    outputbase = hmmscanResults.split(os.sep)[-1].replace(".domtable", "")
    # try to read the domtable file to find out if this gbk has domains. Domains need to be parsed into fastas anyway.
    if os.path.isfile(hmmscanResults):
        pfd_matrix = domtable_parser(outputbase, hmmscanResults)
        num_domains = len(pfd_matrix) # these might still be overlapped, but we need at least 1

        if num_domains > 0:
            print("  Processing domtable file: " + outputbase)

            filtered_matrix, domains = check_overlap(pfd_matrix,overlapCutoff)  #removes overlapping domains, and keeps the highest scoring domain
            
            # Save list of domains per BGC
            pfsoutput = os.path.join(outputdir, outputbase + ".pfs")
            with open(pfsoutput, 'wb') as pfs_handle:
                write_pfs(pfs_handle, domains)
            
            # Save more complete information of each domain per BGC
            pfdoutput = os.path.join(outputdir, outputbase + ".pfd")
            with open(pfdoutput,'wb') as pfd_handle:
                write_pfd(pfd_handle, filtered_matrix)
        else:
            # there aren't any domains in this BGC
            # delete from all data structures
            print("  No domains where found in " + outputbase + ".domtable. Removing it from further analysis")
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
   

def CMD_parser():
    parser = OptionParser()
    
    parser.add_option("-o", "--outputdir", dest="outputdir", default="",
                      help="Output directory, this will contain your pfd, pfs, network and hmmscan output files.")
    parser.add_option("-i", "--inputdir", dest="inputdir", default=os.path.dirname(os.path.realpath(__file__)),
                      help="Input directory of gbk files, if left empty, all gbk files in current and lower directories will be used.")
    parser.add_option("-c", "--cores", dest="cores", default=cpu_count(),
                      help="Set the amount of cores the script may use (default: use all available cores)")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="Prints more detailed information. Toggle to true.")
    parser.add_option("--include_disc_nodes", dest="include_disc_nodes", action="store_true", default=False,
                      help="Include nodes that have no edges to other nodes from the network, default is false. Toggle to true.")
    parser.add_option("-d", "--domain_overlap_cutoff", dest="domain_overlap_cutoff", default=0.1,
                      help="Specify at which overlap percentage domains are considered to overlap")
    parser.add_option("-m", "--min_bgc_size", dest="min_bgc_size", default=0,
                      help="Provide the minimum size of a bgc in base pairs, default is 0bp")
    
    parser.add_option("--seqdist_networks", dest="seqdist_networks", default="A",
                      help="Mode A generates the all vs all networks with sequence distance. Mode S compares clusters within a sample.\
                       Sample input: \"S,A\" generates the samplewise and all vs all. Default is \"A\"")
    
    parser.add_option("--domaindist_networks", dest="domaindist_networks", default="",
                      help="Mode A generates the all vs all networks with domain distance. Mode S compares clusters within a sample.\
                       Sample input: \"A\" only generates the all vs all. Default is \"\"")
    
    parser.add_option("--Jaccardw", dest="Jaccardw", default=0.2,
                      help="Jaccard weight, default is 0.2")
    parser.add_option("--DDSw", dest="DDSw", default=0.75,
                      help="DDS weight, default is 0.75")
    parser.add_option("--GKw", dest="GKw", default=0.05,
                      help="GK weight, default is 0.05")
    parser.add_option("-a", "--anchorboost", dest="anchorweight", default=2.0,
                      help="Boost perceived proportion of anchor DDS subcomponent in 'seqdist' method. Default is to double if (2.0)")
    
    #parser.add_option("--domainsout", dest="domainsout", default="domains",
                      #help="outputfolder of the pfam domain fasta files")
    parser.add_option("--pfam_dir", dest="pfam_dir",
                      default=os.path.dirname(os.path.realpath(__file__)), 
                      help="Location of hmmpress-processed Pfam files. Default is same location of BiG-SCAPE")
    parser.add_option("--anchorfile", dest="anchorfile", default="anchor_domains.txt",
                      help="Provide a custom name for the anchor domains file, default is anchor_domains.txt.")
    parser.add_option("--exclude_gbk_str", dest="exclude_gbk_str", default="",
                      help="If this string occurs in the gbk filename, this file will not be used for the analysis.")
    
    parser.add_option("--mafft_pars", dest="mafft_pars", default="",
                      help="Add single/multiple parameters for mafft specific enclosed by quotation marks e.g. \"--nofft --parttree\"")
    parser.add_option("--al_method", dest="al_method", default="--retree 2",
                      help="alignment method for mafft, if there's a space in the method's name, enclose by quotation marks. default: \"--retree 2\" corresponds to the FFT-NS-2 method")
    parser.add_option("--maxiterate", dest="maxit", default=1000,
                      help="Maxiterate parameter in mafft, default is 1000, corresponds to the FFT-NS-2 method")
    parser.add_option("--mafft_threads", dest="mafft_threads", default=0,
                      help="Set the number of threads in mafft, -1 sets the number of threads as the number of physical cores. Default: same as --cores parameter")
    
    parser.add_option("--force_hmmscan", dest="force_hmmscan", action="store_true", default=False, 
                      help="Force domain prediction using hmmscan even if BiG-SCAPE finds processed domtable files (e.g. to use a new version of PFAM).")
    parser.add_option("--skip_hmmscan", dest="skip_hmmscan", action="store_true", default=False,
                      help="When skipping hmmscan, the GBK files should be available, and the domain tables need to be in the output folder.")
    parser.add_option("--skip_mafft", dest="skip_mafft", action="store_true", default=False, 
                      help="Skip domain prediction by hmmscan as well as domains' sequence's alignments with MAFFT. Needs the original GenBank files, the list of domains per BGC (.pfs) and the BGCs.dict and DMS.dict files.")
    parser.add_option("--skip_all", dest="skip_all", action="store_true",
                      default = False, help = "Only generate new network files. ")
    parser.add_option("--sim_cutoffs", dest="sim_cutoffs", default="1,0.85,0.75,0.6,0.4,0.2",
                      help="Generate networks using multiple similarity (raw distance) cutoff values, example: \"1,0.5,0.1\"")

    parser.add_option("-n", "--nbhood", dest="nbhood", default=4,
                      help="nbhood variable for the GK distance metric, default is set to 4.")

    (options, args) = parser.parse_args()
    return options, args


if __name__=="__main__":
    options, args = CMD_parser()
    
    if options.outputdir == "":
        print "please provide a name for an output folder using parameter -o or --outputdir"
        sys.exit(0)
    
    anchor_domains = get_anchor_domains(options.anchorfile)
    
    global AlignedDomainSequences
    global DomainList
    global verbose
    global BGCs
    global group_dct # contains the class of the Gene Clusters (predicted by antiSMASH and annotated in the GenBank files. Used in the final network files)
    global output_folder
    global pfam_dir
    global timings_file
    global Jaccardw
    global DDSw
    global GKw
    global anchorweight
    global nbhood
    global cores
    include_disc_nodes = options.include_disc_nodes
    cores = int(options.cores)
    nbhood = int(options.nbhood)
    anchorweight = float(options.anchorweight)
    if anchorweight < 1.0:
        sys.exit("Invalid anchorweight parameter (must be equal or greater than 1)")
    Jaccardw = float(options.Jaccardw)
    DDSw = float(options.DDSw) 
    GKw = float(options.GKw)
    cutoff_list = options.sim_cutoffs.split(",")
    if "1" not in cutoff_list:
        cutoff_list.append("1") # compulsory for re-runs
        print("Adding cutoff=1.0 case by default")
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
                    
    verbose = options.verbose
    networks_folder = "networks"
    domainsout = "domains"
    seqdist_networks = options.seqdist_networks.split(",")
    domaindist_networks = options.domaindist_networks.split(",")
    
    if options.skip_all:
        if options.skip_hmmscan or options.skip_mafft:
            print("Overriding --skip_hmmscan/--skip_mafft with --skip_all parameter")
            options.skip_hmmscan = False
            options.skip_mafft = False
    
    time1 = time.time()
    print("\n   - - Obtaining input files - -")
    
    # Make the following available for possibly deleting entries within parseHmmScan
    global genbankDict, gbk_files, sampleDict, clusters, baseNames
    
    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    group_dct = {} # also 
    genbankDict = get_gbk_files(options.inputdir, int(options.min_bgc_size), options.exclude_gbk_str, group_dct)

    # clusters and sampleDict contain the necessary structure for all-vs-all and sample analysis
    clusters = genbankDict.keys()
    
    sampleDict = {} # {sampleName:set(bgc1,bgc2,...)}
    gbk_files = [] # raw list of gbk file locations
    for (cluster, (path, clusterSample)) in genbankDict.iteritems():
        gbk_files.append(path)
        for sample in clusterSample:
            clustersInSample = sampleDict.get(sample, set())
            clustersInSample.add(cluster)
            sampleDict[sample] = clustersInSample

    # This gets very messy if there are files with errors, as the original
    # structures (clusters, sampleDict) are not changed. We'd have to extract
    # only gbk_files, then process it and possibly delete entries, then obtain
    # clusters and sampleDict. genbankDict.keys would still be holding non-valid
    # entries. I'll pass the verifying to get_gbk_files for now
    #print("\nVerifying input files")
    #check_data_integrity([gbk_files]) # turned into a list-within-a-list to retain backwards compatibility
    
    print("\nCreating output directories")
    try:
        os.mkdir(output_folder)
    except OSError as e:
        if "Errno 17" in str(e) or "Error 183" in str(e):
            if not (options.skip_hmmscan or options.skip_all or options.skip_mafft):
                print(" Warning: Output directory already exists!")
            else:
                print(" Using existing output directory.")
        else:
            print("Unexpected error when creating output directory")
            sys.exit(str(e))
    write_parameters(output_folder, options)
    
    try:
        os.mkdir(os.path.join(output_folder, networks_folder))
    except OSError as e:
        if "Errno 17" in str(e) or "Error 183" in str(e):
            print(" Warning: possibly overwriting files in network folder")
            pass
        else:
            print("Unexpected error when creating network directory")
            sys.exit(str(e))
    
    try:
        os.mkdir(os.path.join(output_folder, domainsout))
    except OSError as e:
        # 17 (Linux): "[Errno 17] File exists";
        # 183 (Windows) "[Error 183] Cannot create a file when that file already exists"
        if "Errno 17" in str(e) or "Error 183" in str(e):
            if not (options.skip_all or options.skip_hmmscan or options.skip_mafft):
                print(" Emptying domains directory")
                for thing in os.listdir(os.path.join(output_folder, domainsout)):
                    os.remove(os.path.join(output_folder, domainsout, thing))
            else:
                print(" Using existing domains directory")
        else:
            print("Fatal error when trying to create domains' directory")
            sys.exit(str(e))


    if verbose:
        print(" Trying threading on %i cores" % cores)
    
    #open the file that will contain the timed functions
    timings_file = open(os.path.join(output_folder, "runtimes.txt"), 'w') 
    
    """BGCs -- 
    dictionary of this structure:
    BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""     
    BGCs = {} #will contain the BGCs
    
    AlignedDomainSequences = {} # Key: specific domain sequence label. Item: aligned sequence
    DomainList = {} # Key: BGC. Item: ordered list of domains
    
    # to avoid calling MAFFT if there's only 1 seq. representing a particular domain
    sequences_per_domain = {} 
    
    print("\n\n   - - Processing input files - -")
    
    # These will be used to track if we drop files in processing
    genbankFileLocations = set(gbk_files)
    baseNames = set(clusters)

    ### Step 1: Generate Fasta Files
    print "\nParsing genbank files to generate fasta files for hmmscan"

    # filter through task list to avoid unecessary computation: if output file is already there and non-empty, exclude it from list
    alreadyDone = set()
    for genbank in genbankFileLocations:
        outputbase = genbank.split(os.sep)[-1].replace(".gbk","")
        outputfile = os.path.join(output_folder,outputbase + '.fasta')
        if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
            alreadyDone.add(genbank)

    if len(genbankFileLocations - alreadyDone) == 0:
        print(" All GenBank files had already been processed")
    elif len(alreadyDone) > 0:
        print " Warning! The following NEW input file(s) will be processed: %s" % ", ".join(x.split(os.sep)[-1].replace(".gbk","") for x in genbankFileLocations - alreadyDone)
    else:
        print(" Processing " + str(len(genbankFileLocations)) + " files")

    # Generate Pool of workers
    pool = Pool(cores,maxtasksperchild=32)
    for genbankFile in (genbankFileLocations - alreadyDone):
        pool.apply_async(generateFasta,args =(genbankFile,output_folder))
    pool.close()
    pool.join()
    print " Finished generating fasta files."


    ### Step 2: Run hmmscan
    print("\nPredicting domains using hmmscan")
    
    # All available fasta files (could be more than it should if reusing output folder)
    allFastaFiles = set(glob(os.path.join(output_folder,"*.fasta")))
    
    # fastaFiles: all the fasta files that should be there 
    # (i.e. correspond to the input files)
    fastaFiles = set()
    for name in baseNames:
        fastaFiles.add(os.path.join(output_folder, name+".fasta"))
    
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
            outputbase  = fasta.split(os.sep)[-1].replace(".fasta","")
            outputfile = os.path.join(output_folder,outputbase + '.domtable')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(fasta)
            
        if len(fastaFiles - alreadyDone) == 0:
                print(" All fasta files had already been processed")
        elif len(alreadyDone) > 0:
            print " Warning! The following NEW fasta file(s) will be processed: %s" % ", ".join(x.split(os.sep)[-1].replace(".fasta","") for x in fastaFiles - alreadyDone)
        else:
            print(" Predicting domains for " + str(len(fastaFiles)) + " fasta files")

        task_set = fastaFiles - alreadyDone
        
    pool = Pool(cores,maxtasksperchild=1)
    for fastaFile in task_set:
        pool.apply_async(runHmmScan,args=(fastaFile,pfam_dir,output_folder, verbose))
    pool.close()
    pool.join()

    print " Finished generating domtable files."


    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    print("\nParsing hmmscan domtable files")
    
    # All available domtable files
    allDomtableFiles = set(glob(os.path.join(output_folder,"*.domtable")))
    
    # domtableFiles: all domtable files corresponding to the input files
    domtableFiles = set()
    for name in baseNames:
        domtableFiles.add(os.path.join(output_folder, name+".domtable"))
    
    # domtableBases: the actual set of input files with coresponding domtable files
    domtableBases = allDomtableFiles.intersection(domtableFiles)
    
    # Verify that all input files have a corresponding domtable file
    if len(domtableFiles - domtableBases) > 0:
        sys.exit("Error! The following files did NOT have their domains predicted: " + ", ".join(domtableFiles - domtableBases))
    
    # find already processed files
    alreadyDone = set()
    if not options.force_hmmscan:
        for domtable in domtableFiles:
            outputbase = domtable.split(os.sep)[-1].replace(".domtable","")
            outputfile = os.path.join(output_folder,outputbase + '.pfd')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(domtable)

    if len(domtableFiles - alreadyDone) == 0:
        print(" All domtable files had already been processed")
    elif len(alreadyDone) > 0:
        print " Warning! The following domtable files had not been processed: %s" % ", ".join(x.split(os.sep)[-1].replace(".domtable","") for x in domtableFiles - alreadyDone)
    else:
        print(" Processing " + str(len(domtableFiles)) + " domtable files")

    # If using the multiprocessing version and outputbase doesn't have any
    #  predicted domains, it's not as easy to remove if from the analysis
    #  (probably because parseHmmScan only has a copy of clusters et al?)
    # Using serialized version for now. Probably doesn't have too bad an impact
    #pool = Pool(cores,maxtasksperchild=32)
    for domtableFile in domtableFiles - alreadyDone:
        parseHmmScan(domtableFile,output_folder,options.domain_overlap_cutoff)
        #pool.apply_async(parseHmmScan, args=(domtableFile,output_folder,options.domain_overlap_cutoff))
    #pool.close()
    #pool.join()

    print " Finished generating generating pfs and pfd files."


    ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per sample objects
    print("\nProcessing domains sequence files")
    
    # All available pfd files
    allPfdFiles = set(glob(os.path.join(output_folder,"*.pfd")))
    
    # pfdFiles: all pfd files corresponding to the input files
    # (some input files could've been removed due to not having predicted domains)
    pfdFiles = set()
    for name in baseNames:
        pfdFiles.add(os.path.join(output_folder, name+".pfd"))
    
    # pfdBases: the actual set of input files that have pfd files
    pfdBases = allPfdFiles.intersection(pfdFiles)
    
    # verify previous step. 
    # All BGCs without predicted domains should no longer be in baseNames    
    if len(pfdFiles - pfdBases) > 0:
        sys.exit("Error! The following files did NOT have their domtable files processed: " + ", ".join(pfdFiles - pfdBases))

    if options.skip_mafft:
        print(" Running with skip_mafft parameter: Assuming that the domains folder has all the fasta files")
        print(" Only extracting BGC group from input file")
    else:
        if verbose:
            print(" Adding sequences to corresponding domains file")
            
        for outputbase in baseNames:
            if verbose:
                print("   Processing: " + outputbase)

            pfdFile = os.path.join(output_folder, outputbase + ".pfd")
            filtered_matrix = [map(lambda x: x.strip(), line.split('\t')) for line in open(pfdFile)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(output_folder, outputbase + ".fasta")
            fasta_dict = fasta_parser(open(fasta_file, "r")) # all fasta info from a BGC
            save_domain_seqs(filtered_matrix, fasta_dict, domainsout, output_folder, outputbase)

            BGCs[outputbase] = BGC_dic_gen(filtered_matrix)

    # Get the ordered list of domains
    print(" Reading the ordered list of domains from the pfs files")
    for outputbase in baseNames:
        pfsfile = os.path.join(output_folder, outputbase + ".pfs")
        if os.path.isfile(pfsfile):
            DomainList[outputbase] = get_domain_list(pfsfile)
        else:
            sys.exit(" Error: could not open " + outputbase + ".pfs")

    #Write or retrieve BGC dictionary
    if not options.skip_all:
        if options.skip_hmmscan or options.skip_mafft:
            with open(os.path.join(output_folder, "BGCs.dict"), "r") as BGC_file:
                BGCs = pickle.load(BGC_file)
                BGC_file.close()
        else:
            with open(os.path.join(output_folder, "BGCs.dict"), "w") as BGC_file:
                pickle.dump(BGCs, BGC_file)
                BGC_file.close()
    
    
    print("\n\n   - - Calculating distance matrix - -")
    
    
    # Distance without taking sequence similarity between specific domains into account
    if domaindist_networks:
        if options.skip_all: #read already calculated distances
            print(" Trying to read alread calculated network file...")
            if os.path.isfile(os.path.join(output_folder, networks_folder, "networkfile_domain_dist_all_vs_all_c1.network")):
                network_matrix = network_parser(os.path.join(output_folder, networks_folder, "networkfile_domain_dist_all_vs_all_c1.network"), Jaccardw, DDSw, GKw, anchorweight)
                print("  ...done")
            else:
                sys.exit("  File networkfile_domain_dist_all_vs_all_c1.network could not be found!")
            
        if 'A' in domaindist_networks:
            print("\nGenerating all-vs-all network with domain distance method")
            pairs = set(map(tuple, map(sorted, combinations(clusters, 2))))
            cluster_pairs = [(x, y, "domain_dist", anchor_domains) for (x, y) in pairs]
            network_matrix = generate_network(cluster_pairs, cores)
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, os.path.join(output_folder, networks_folder, "networkfile_domain_dist_all_vs_all_c" + cutoff + ".network"), include_disc_nodes)
            if 'S' in domaindist_networks:
                if len(sampleDict) == 1:
                    print("\nNOT generating networks per sample (only one sample, covered in the all-vs-all case)")
                else:
                    print("\nGenerating sample networks with domain distance method")
                    for sample, sampleClusters in sampleDict.iteritems():
                        print(" Sample: " + sample)
                        if len(sampleClusters) == 1:
                            print(" Warning: Sample size = 1 detected. Not generating networks for this sample (" +
                                  sample + ")")
                        else:
                            pairs = set(map(tuple, map(sorted, combinations(sampleClusters, 2))))
                            network_matrix_sample = {}
                            for pair in pairs:
                                network_matrix_sample[pair] = network_matrix[pair]
                            for cutoff in cutoff_list:
                                write_network_matrix(network_matrix_sample, cutoff,
                                                     os.path.join(output_folder, networks_folder,
                                                                  "networkfile_domain_dist_" + sample + "_c" + cutoff + ".network"),
                                                     include_disc_nodes)
        elif 'S' in domaindist_networks:
            # need to caculate the network for each of the pairs
            if len(sampleDict) == 1:
                print("\nNOT generating networks per sample (only one sample, covered in the all-vs-all case)")
            else:
                print("\nGenerating sample networks with domain distance method")
                for sample, sampleClusters in sampleDict.iteritems():
                    print(" Sample: " + sample)
                    if len(clusters) == 1:
                        print(" Warning: Sample size = 1 detected. Not generating networks for this sample (" +
                              sample + ")")
                    else:
                        pairs = set(map(tuple, map(sorted, combinations(sampleClusters, 2))))
                        cluster_pairs = [(x, y, "domain_dist", anchor_domains) for (x, y) in pairs]
                        network_matrix_sample = generate_network(cluster_pairs, cores)
                        for cutoff in cutoff_list:
                            write_network_matrix(network_matrix_sample, cutoff,
                                                 os.path.join(output_folder, networks_folder,
                                                              "networkfile_domain_dist_" + sample + "_c" + cutoff + ".network"),
                                                 include_disc_nodes)
                            # Need to calculate the networks per sample from the all-v-all network matrix
    # Check whether user wants seqdist method networks before calculating DMS

    if seqdist_networks:
        if options.skip_all:
            print(" Trying to read already calculated network file...")
            if os.path.isfile(os.path.join(output_folder, networks_folder, "networkfile_seqdist_all_vs_all_c1.network")):
                network_matrix = network_parser(os.path.join(output_folder, networks_folder, "networkfile_seqdist_all_vs_all_c1.network"), Jaccardw, DDSw, GKw, anchorweight)
                print("  ...done")
            else:
                sys.exit("  File networkfile_seqdist_all_vs_all_c1.network could not be found!")
            
        elif not options.skip_mafft:
            # obtain all fasta files with domain sequences
            fasta_domains = get_domain_fastas(domainsout, output_folder)

            for domain_file in fasta_domains:
                domain_name = domain_file.split(os.sep)[-1].replace(".fasta", "")
                
                # fill fasta_dict...
                with open(domain_file, "r") as fasta_handle:
                    fasta_dict = fasta_parser(fasta_handle)
                # ...to find out how many sequences do we actually have
                if len(fasta_dict) == 1:
                    # avoid calling MAFFT if it's not possible to align (only one sequence)
                    if verbose:
                        print(" Skipping MAFFT for domain " + domain_name + " (only one sequence)")
                else:           
                    if verbose:
                        print(" Running MAFFT for domain: " + domain_name)
                    
                    domain_file_base = domain_file.replace(".fasta", "")
                    
                    # Multiple alignment of all domain sequences
                    run_mafft(options.al_method, options.maxit, options.mafft_threads, options.mafft_pars, domain_file_base)
                    
                    # Check if MAFFT's output file was generated
                    if not os.path.isfile(domain_file_base + ".algn"):
                        print("  Warning, " + domain_name + ".algn could not be found (did MAFFT failed?)")
                    
  
        print(" Trying to read domain alignments (*.algn files)")            
        aligned_files_list = glob(os.path.join(output_folder, domainsout, "*.algn"))
        if len(aligned_files_list) == 0:
            sys.exit("No aligned sequences found in the domain folder (run without the --skip_mafft parameter or point to the correct output folder)")
        for aligned_file in aligned_files_list:
            with open(aligned_file, "r") as aligned_file_handle:
                fasta_dict = fasta_parser(aligned_file_handle)
                for header in fasta_dict:
                    AlignedDomainSequences[header] = fasta_dict[header]
            
        if "A" in seqdist_networks:
            print("\nGenerating all-vs-all network with domain-sequence distance method")
            if not options.skip_all:
                print(" Calculating all pairwise distances")
                pairs = set(map(tuple, map(sorted, combinations(clusters, 2))))
                cluster_pairs = [(x, y, "seqdist", anchor_domains) for (x, y) in pairs]
                network_matrix = generate_network(cluster_pairs, cores)
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, os.path.join(output_folder, networks_folder, "networkfile_seqdist_all_vs_all_c" + cutoff + ".network"), include_disc_nodes)
                
        if "S" in seqdist_networks:
            if len(sampleDict) == 1 and "A" in seqdist_networks:
                print("\nNOT generating networks per sample (only one sample, covered in the all-vs-all case)")
            else:
                print("\nGenerating sample networks with domain-sequence distance method")
                for sample, sampleClusters in sampleDict.iteritems():
                    print(" Sample: " + sample)
                    if len(sampleClusters) == 1:
                        print(" Warning: Sample size = 1 detected. Not generating network for this sample (" + sample + ")")
                    else:
                        pairs = set(map(tuple, map(sorted, combinations(sampleClusters, 2))))
                        cluster_pairs = [(x, y, "seqdist", anchor_domains) for (x, y) in pairs]
                        network_matrix_sample = {}
                        if "A" in seqdist_networks or options.skip_all:
                            for pair in pairs:
                                network_matrix_sample[pair] = network_matrix[pair]
                        else:
                            network_matrix_sample = generate_network(cluster_pairs, cores)
                        for cutoff in cutoff_list:
                            write_network_matrix(network_matrix_sample, cutoff,
                                                 os.path.join(output_folder, networks_folder,
                                                              "networkfile_seqdist_" + sample + "_c" + cutoff + ".network"),
                                                 include_disc_nodes)

    runtime = time.time()-time1
    runtime_string = '\tMain function took %0.3f s' % (runtime)
    timings_file.write(runtime_string + "\n")
    print runtime_string
    
    timings_file.close()
