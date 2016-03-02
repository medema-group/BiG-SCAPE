#!/usr/bin/env python

"""
Student/programmer: Marley Yeong
marleyyeong@live.nl
supervisor: Marnix Medema
dependencies: hmmer, biopython, mafft, munkres, numpy

Usage:   place this script in the same folder as the antismash output files
         include the munkres.py file in this folder
         include the files necessary for hmmscan in this folder
         $python bgc_networks.py

Status: development/testing


Todo: 
implement anchor domains
"""

import fileinput, pickle, sys, math
from munkres import Munkres
import numpy as np
import sys
import re
import os
import subprocess
import time
from optparse import OptionParser
import math
import urllib2 
import multiprocessing
from functools import partial

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from site import setBEGINLIBPATH


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

def timeit(f):
    def wrap(*args):
        insignificant_runtime = 1
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        runtime = time2-time1
        
        runtime_string = '%s function took %0.3f s' % (f.func_name, runtime)
        #write_timings(runtime_string)
        #=======================================================================
        # timings_file.write(runtime_string + "\n")
        #=======================================================================
        if runtime > insignificant_runtime:
            
            print runtime_string
            
        return ret
    return wrap


@timeit
def distout_parser(distout_file):
    """returns distance values, for domains in the following format  { ('specific_domain_name_1',
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
                
                #return {('', ''): (0.000000001, 0)}
            
            seqsdict[seq_number] = "".join(line.split("=")[1:]).strip()#in case the header contains an = sign

        elif linecounter > 3 + numberof_seqs:
            distances += line.strip().split(" ")

    keys=[]
    if len(distances) != (numberof_seqs * (numberof_seqs-1)) / 2.0:
        print "something went horribly wrong in importing the distance matrix"
    else:
        print "distance matrix imported correctly"
        keys = seqsdict.keys()
        print sorted(keys)

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


@timeit
def iterFlatten(root):
    if isinstance(root, (list, tuple)):
        for element in root:
            for e in iterFlatten(element):
                yield e
    else:
        yield root


def write_network_matrix(matrix, cutoff, filename, exclude_disc_nodes):
    networkfile = open(filename, 'w+')
    clusters = [] # will contain the names of clusters that have an edge value lower than the threshold
    networkfile.write("clustername1\tclustername2\tgroup1\tdefinition\tgroup2\tdefinition\t-log2score\traw distance\tsquared similarity\tcombined group\tshared group\n")
    for row in matrix:
        temprow = []
        #if both clusters have the same group, this group will be annotated in the last column of the network file
        #print "raw_distance", float(row[5]) 
        #print "cutoff", float(cutoff)
        for i in row:
            temprow.append(i)
        
        if float(temprow[7]) <= float(cutoff):
            clusters.append(row[0])
            clusters.append(row[1])

            if row[2] != "" and row[3] != "":
                temprow.append(" - ".join(sorted([str(row[2]),str(row[3])])))
            elif row[3] != "":
                temprow.append(str(row[3]))
            elif row[2] != "":
                temprow.append(str(row[2]))
            else:
                 temprow.append(str("NA"))
            
            if row[2] == row[3]:
                temprow.append(row[2])
            else:
                temprow.append("")
                
            networkfile.write("\t".join(temprow) + "\n")

    if exclude_disc_nodes == False:   
        clusters = set(clusters)
    #clustername1    clustername2    group1    definition    group2    definition    #NAME?    raw distance    squared similarity    combined group    shared group
        passed_clusters = []
        for row in matrix:
            if row[0] not in clusters and row[0] not in passed_clusters:
                networkfile.write("\t".join([row[0],row[0],row[2],row[3],'','',str(0),str(0),str(0),'','']) + "\n")
                passed_clusters.append(row[0])
            elif row[1] not in clusters and row[1] not in passed_clusters:
                networkfile.write("\t".join([row[1],row[1],row[4],row[5],'','',str(0),str(0),str(0),'','']) + "\n")
                passed_clusters.append(row[1])
            
    networkfile.close()   



def generate_dist_matrix(parms):
    
    #Get the values from the parameters
    cluster1 = parms[0]
    cluster2 = parms[1]
    dist_method = parms[2]
    anchor_domains = parms[3]
    
    cluster_file1 = outputdir + "/" + cluster1 + ".pfs"
    cluster_file2 = outputdir + "/" + cluster2 + ".pfs"
    if dist_method == "seqdist":
        dist = cluster_distance(cluster1, cluster2, anchor_domains)
    elif dist_method == "domain_dist":
        dist = Distance_modified(get_domain_list(cluster_file1), get_domain_list(cluster_file2), 0, 4)
        
    if dist == 0:
        logscore = float("inf")
    else:
        logscore = 0
        try:
            logscore = -math.log(dist, 2) #Write exception, ValueError
        except ValueError:
            print "calculating the logscore with distance", dist, "failed in function generate_network, using distance method", dist_method, networkfilename
            
        
    #clustername1 clustername2 group1, group2, -log2score, dist, squared similarity
    network_row = [str(cluster1), str(cluster2),group_dct[cluster1][0],group_dct[cluster1][1],\
    group_dct[cluster2][0],group_dct[cluster2][1], str(logscore), str(dist), str((1-dist)**2)]
    
    return network_row
    
        
@timeit
def generate_network(bgc_list, group_dct, networkfilename, networks_folder, dist_method, anchor_domains):
    #Contents of the network file: clustername1 clustername2 group1, group2, -log2score, dist, squared similarity
    "saves the distances as the log2 of the similarity"
    network_matrix_conv = []

    
    if networkfilename == "":
    
        networkfilename = outputdir + "/" + networks_folder + "/" + "networkfile_" + dist_method + "_" + \
        str("".join(bgc_list[0].split(".")[0:-2]))
    else:
        networkfilename = outputdir + "/" + networks_folder + "/" + "networkfile_" + dist_method + "_" + \
        networkfilename
        
    print "Generating network:", networkfilename
    
    print "*************", bgc_list

    cluster_queue = []
    cluster_pairs = []
    
    for bgc in bgc_list:
        cluster_queue.append(bgc) #deep copy
        
    for cluster1 in bgc_list:
        cluster_queue.remove(cluster1)
        
        for cluster2 in cluster_queue:
            #addclusterpair
            cluster_pairs.append([cluster1, cluster2, dist_method, anchor_domains])    
            
    
    
    pool = multiprocessing.Pool(cores) #create the appropriate amount of pool instances
    network_matrix = pool.map(generate_dist_matrix, cluster_pairs)  #Assigns the data to the different workers and pools the results
                                                                    #back into the network_matrix variable   

    sorted_network_matrix = sorted(network_matrix, key=lambda network_matrix_entry: network_matrix_entry[5]) #sort the matrix on the log2 scores
    return sorted_network_matrix, networkfilename
    

def cluster_distance(A,B, anchor_domains): 
  
    #key is name of GC, values is list of specific pfam domain names
    clusterA = BGCs[A] #will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    clusterB = BGCs[B] #{ 'general_domain_name_x' : ['specific_domain_name_1', 'specific_domain_name_2'], etc }
    #A and B are lists of pfam domains
     
    try:
        #calculates the intersect
        Jaccard = len(set(clusterA.keys()) & set(clusterB.keys())) / \
              float( len(set(clusterA.keys())) + len(set(clusterB.keys())) \
              - len(set(clusterA.keys()) & set(clusterB.keys())))
    except ZeroDivisionError:
        print "Zerodivisionerror during the Jaccard distance calculation. Can only happen when both clusters are empty"
        print "keys of clusterA", A, clusterA.keys()
        print "keys of clusterB", B, clusterB.keys()
         
    intersect = set(clusterA.keys() ).intersection(clusterB.keys()) #shared pfam domains
    not_intersect = []
    #DDS: The difference in abundance of the domains per cluster
    #S: Max occurence of each domain
    DDSa,Sa = 0,0
    DDS,S = 0,0
    SumDistance = 0
    pair = ""

    for domain in set(clusterA.keys() + clusterB.keys()):
        if domain not in intersect:
            not_intersect.append(domain)
            
    for unshared_domain in not_intersect: #no need to look at seq identity or anchors, since these domains are unshared
        #for each occurence of an unshared domain do DDS += count of domain and S += count of domain
        dom_set = []
        try:
            dom_set = clusterA[unshared_domain]
        except KeyError:
            dom_set = clusterB[unshared_domain]
            
        DDS += len(dom_set)
        S += len(dom_set)
        
        
    for shared_domain in intersect:
        seta = clusterA[shared_domain]
        setb = clusterB[shared_domain]
        
        if len(seta+setb) == 2: #The domain occurs only once in both clusters
            pair = tuple(sorted([seta[0],setb[0]]))
            SumDistance = 1-DMS[shared_domain][pair][0] #1
            if shared_domain.split(".")[0] in anchor_domains: 
                Sa += max(len(seta),len(setb))
                DDSa += SumDistance 
            else:
                S += max(len(seta),len(setb))
                DDS += SumDistance
        else:                   #The domain occurs more than once in both clusters
            accumulated_distance = 0
            
            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            for domsa in range(len(seta)):
                for domsb in range(domsa, len(setb)):
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    Similarity = DMS[shared_domain][pair][0]
                    
                    seq_dist = 1-Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist
                
            #Only use the best scoring pairs
            Hungarian = Munkres()
            #print "DistanceMatrix", DistanceMatrix
            BestIndexes = Hungarian.compute(DistanceMatrix)
            #print "BestIndexes", BestIndexes
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            #print "accumulated_distance", accumulated_distance
                         
                    
            SumDistance = (abs(len(seta)-len(setb)) + accumulated_distance)  #diff in abundance + sequence distance
            
            if shared_domain.split(".")[0] in anchor_domains: 
                Sa += max(len(seta),len(setb))
                DDSa += SumDistance 
            else:
                S += max(len(seta),len(setb))
                DDS += SumDistance

 
    #  calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A]
    Ar.reverse()
    GK = max([calculate_GK(A, B, nbhood), calculate_GK(Ar, B, nbhood)])
     
    if DDSa != 0:
        DDSa /= float(Sa)
        DDS = (anchorweight * DDSa) + (1 - anchorweight) * DDS    #Recalculate DDS by giving preference to anchor domains
        DDS /= float(S)
    else:
        DDS /= float(S) 
    
    
    DDS = 1-DDS #transform into similarity
    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK) 
    if Distance < 0:
        print "negative distance", Distance, "DDS", DDS, pair
        print "Probably a rounding issue"
        print "Distance is set to 0 for these clusters" 
        Distance = 0
    return Distance



# --- compare two clusters using distance information between sequences of domains
@timeit
def CompareTwoClusters(A,B): 
    anchor_domains = ["PF00668"]
  
    #key is name of GC, values is list of specific pfam domain names
    clusterA = BGCs[A] #will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    clusterB = BGCs[B] #{ 'general_domain_name_x' : ['specific_domain_name_1', 'specific_domain_name_2'], etc }
    #A and B are lists of pfam domains
     
    try:
        #calculates the intersect
        Jaccard = len(set(clusterA.keys()) & set(clusterB.keys())) / \
              float( len(set(clusterA.keys())) + len(set(clusterB.keys())) \
              - len(set(clusterA.keys()) & set(clusterB.keys())))
    except ZeroDivisionError:
        print "Zerodivisionerror during the Jaccard distance calculation. Can only happen when both clusters are empty"
        print "keys of clusterA", A, clusterA.keys()
        print "keys of clusterB", B, clusterB.keys()
         
    #DDS: The difference in abundance of the domains per cluster
    #S: Max occurence of each domain
    DDSa,Sa = 0,0
    DDS,S = 0,0
    anchor = False
    SumDistance = 0
    for domain in set(clusterA.keys() + clusterB.keys()):
         
        try: seta = clusterA[domain] #Get the specific domains, that can be mapped to the sequence identity
        except KeyError: 
            seta = []
             
        try: setb = clusterB[domain]
        except KeyError: 
            setb = []
            


        if len(seta) == 1 and len(setb) == 1: #if both clusters only contain this domain once
            pair = tuple(sorted([seta[0],setb[0]]))
            try: SumDistance = 1-DMS[domain][pair][0] #1
            except KeyError: 
                print "Could not retrieve the distance metric for this cluster pair", pair
                SumDistance = 1-0.1
          #  print "if", SumDistance
            if domain.split(".")[0] in anchor_domains:
                
                Sa += max(len(seta),len(setb))
                DDSa += SumDistance #Also include similar code for the if/elif statements above
            else:
                S += max(len(seta),len(setb))
                DDS += SumDistance
 
        elif len(seta) + len(setb)==1: #If one cluster does not contain this domain
            DDS += 1.
            S += 1.
 
        else:                           #if the domain occurs more than once in either or both of the clusters
            print "seta", seta
            print "setb", setb
            N = max(len(seta),len(setb))
            DistanceMatrix = np.zeros((N,N)) #creates a matrix of N by N dimensions filled with zeros
            for ja in xrange(len(seta)):
                for jb in xrange(ja,len(setb)):
                    pair = tuple(sorted([seta[ja],setb[jb]]))
                    try: 
                        Distance = 1-DMS[domain][pair][0]
                    except KeyError: 
                        Distance = 1-0.1
                        print "KeyError when retrieving the distance score during the hungarian calculation", pair
                      
                    print "Distance", Distance
                      
                    DistanceMatrix[ja,jb] = Distance
                    DistanceMatrix[jb,ja] = Distance
 
 
            print "DistanceMatrix", DistanceMatrix
            Hungarian = Munkres()
            BestIndexes = Hungarian.compute(DistanceMatrix)
            SumDistance = sum([DistanceMatrix[bi] for bi in BestIndexes])
 
            if domain.split(".")[0] in anchor_domains:
                Sa += max(len(seta),len(setb))
                DDSa += SumDistance
            else:
                S += max(len(seta),len(setb))
                DDS += SumDistance
 
    #  calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A]
    Ar.reverse()
    GK = max([calculate_GK(A, B, nbhood), calculate_GK(Ar, B, nbhood)])
     
    if DDSa != 0:
        DDSa /= float(Sa)
        DDS = (anchorweight * DDSa) + (1 - anchorweight) * DDS    #Recalculate DDS by giving preference to anchor domains
        #print "DDSa", DDS
        DDS /= float(S)
        #print "DDS", DDS
    else:
        DDS /= float(S) 
     
    DDS = math.exp(-DDS) #will 'flip' the DDS value 0.9 becomes 0.4, 0.1 becomes 0.9 
     
    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK) 
   # Distance = 1 - 0.36*Jaccard - 0.64*DDS 0.63
         
    #===========================================================================
    # Similarity = 1-Distance
    # lin = '%s\t%s\t%.4f\n' % (A,B,Similarity)
    # output = open("seqdist.txt", 'w')
    # output.write(lin)
    #===========================================================================
    
    return Distance


@timeit
def run_mafft(al_method, maxit, threads, mafft_pars, domain):
    """Runs mafft program with the provided parameters.
    The specific parameters used to run mafft with are actually very important for the final result.
    Using mafft with the most progressive parameters does indeed affect the quality of the distance.
    It is better to just use the domain information for the distance if more computationally intensive options
    for mafft cannot be used. Setting maxiterate to anything higher than 10 did not make a difference in accuracy in the nrps testset"""
    
    alignment_file = domain + ".algn"
    
    
    mafft_cmd_list = []
    mafft_cmd_list.append("mafft --distout")
    mafft_cmd_list.append("--quiet")
    mafft_cmd_list.append(al_method)
    if maxit != 0:
        mafft_cmd_list.append("--maxiterate " + str(maxit))
        
    mafft_cmd_list.append("--thread")
    mafft_cmd_list.append(str(threads))
    
    if mafft_pars != "":
        mafft_cmd_list.append(mafft_pars)
        
    mafft_cmd_list.append(str(domain) + ".fasta")
    mafft_cmd_list.append(">")
    mafft_cmd_list.append("alignment_file")
    
    mafft_cmd = " ".join(mafft_cmd_list)
    
    print mafft_cmd
    subprocess.check_output(mafft_cmd, shell=True)


@timeit
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

@timeit
def get_domain_fastas(domain_folder):
    """Finds the pfam domain fasta files"""
    domain_fastas = []
    for dirpath, dirnames, filenames in os.walk(outputdir + "/" + domain_folder + "/"):
        for fname in filenames:
            if ".fasta" in fname and "hat2" not in fname:
                domain_fastas.append(dirpath + "/" + fname)
                if verbose == True:
                    print fname
                    
    return domain_fastas

@timeit    
def calc_perc_identity(seq1, seq2):
    """Percent Identity = (Matches x 100)/Length of aligned region (with gaps)
    Note that only internal gaps are included in the length, not gaps at the sequence ends."""

    al_len = min([len(seq1.strip("-")), len(seq2.strip("-"))])
    matches = 0
    for pos in range(len(seq1)): #Sequences should have the same length because they come from an MSA
        try:
            if seq1[pos] == seq2[pos] and not (seq1[pos] == "-" and seq2[pos] == "-"):  #don't count two gap signs as a match!
                matches += 1
        except IndexError:
            print "Something went wrong, most likely your alignment file contains duplicate sequences"
            

    return (matches * 100) / float(al_len), al_len    


@timeit
def BGC_dic_gen(filtered_matrix):
    """Generates the: { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } part of the BGCs variable."""
    bgc_dict = {}
    for row in filtered_matrix:
        try: #Should be faster than performing if key in dictionary.keys()
            bgc_dict[row[6]]
            bgc_dict[row[6]].append(str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4]))
            #bgc_dict[row[6]].append(str(row[6]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4]))
        except KeyError: #In case of this error, this is the first occurrence of this domain in the cluster
           bgc_dict[row[6]]=[str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4])]
            #bgc_dict[row[6]]=[str(row[6]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4])]
            
    return bgc_dict

    
@timeit
def save_domain_seqs(filtered_matrix, fasta_dict, domains_folder):
    """Write fasta sequences for the domains in the right pfam-domain file"""
    for row in filtered_matrix:
        domain = row[6]
        seq = fasta_dict[">"+str(row[-1].strip())] #access the sequence by using the header
        domain_file = open(outputdir + "/" + domains_folder + "/" + domain +".fasta", 'a')
        #same as in BGCs variable
        domain_file.write(">" + str(row[0]) + "_" + str(row[-1]) + "_" + str(row[3]) + "_" + str(row[4]) \
        + "\n" + str(seq)[int(row[3]):int(row[4])] + "\n")
        domain_file.close()
        

@timeit
def calculate_GK(A, B, nbhood):
    # calculate the Goodman-Kruskal gamma index
    GK = 0.
    if len(set(A) & set(B)) > 1:
        pairsA = set( [(A[i],A[j]) for i in xrange(len(A)-nbhood) for j in xrange(i+1,i+nbhood)] )
        pairsB = set( [(B[i],B[j]) for i in xrange(len(B)-nbhood) for j in xrange(i+1,i+nbhood)] )
        allPairs = set(list(pairsA) + list(pairsB))
        Ns, Nr = 0.,0.
        for p in allPairs:
            if p in pairsA and p in pairsB: Ns += 1
            elif p in pairsA and tuple(p[::-1]) in pairsB: Nr += 1
            elif tuple(p[::-1]) in pairsA and p in pairsB: Nr += 1
            else: pass
        if (Nr + Ns) == 0:
            gamma = 0
        else:
            gamma = abs(Nr-Ns) / (Nr+Ns)
        GK = (1+gamma)/2.
    return GK


@timeit
def Distance_modified(clusterA, clusterB, repeat=0, nbhood=4):
    "Modified to work better for 'OBU' detection"
    "Kui Lin, Lei Zhu and Da-Yong Zhang "

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
    DDS = math.exp(-DDS) #transforms the DDS to a value between 0 - 1

    # calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A]
    Ar.reverse()
    GK = max([calculate_GK(A, B, nbhood), calculate_GK(Ar, B, nbhood)])


    # calculate the distance
    #print "Jaccard", Jaccard
    #print "DDS", DDS
    #print "GK", GK
    Distance = 1 - Jaccardw*Jaccard - DDSw*DDS - GKw*GK
    
    if Distance < 0:
        Distance = 0

    return Distance

@timeit
def no_overlap(locA1, locA2, locB1, locB2):
    """Return True if there is no overlap between two regions"""
    if locA1 < locB1 and locA2 < locB1:
        return True
    elif locA1 > locB2 and locA2 > locB2:
        return True
    else:
        return False

@timeit
def overlap_perc(overlap, len_seq):
    return float(overlap) / len_seq
    
    
@timeit
def overlap(locA1, locA2, locB1, locB2):

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


@timeit
def check_overlap(pfd_matrix, overlap_cutoff):
    """Check if domains overlap for a certain overlap_cutoff.
     If so, remove the domain(s) with the lower score."""
     
##pfd example row:
##AB050629.gbk	score   yxjC	1	1167	+	PF03600	CitMHS

    #pfd_matrix = sorted(pfd_matrix, key=lambda pfd_matrix_entry: pfd_matrix_entry[3]) #sort the domains in the cluster by their first coordinate in the 4th col

    row1_count = 0
    delete_list = []
    for row1 in pfd_matrix:
        row1_count += 1
        row2_count = 0
        for row2 in pfd_matrix:
            row2_count += 1
            if row1_count != row2_count and row1[8] == row2[8]: #check if we are not comparing the same rows, but has the same CDS
                #check if there is overlap between the domains
                if no_overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4])) == False:
                    overlapping_nucleotides = overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4]))
                    overlap_perc_loc1 = overlap_perc(overlapping_nucleotides, int(row1[4])-int(row1[3]))
                    overlap_perc_loc2 = overlap_perc(overlapping_nucleotides, int(row2[4])-int(row2[3]))
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[1]) > float(row2[1]):
                            delete_list.append(row2)
                        elif float(row1[1]) < float(row2[1]):
                             delete_list.append(row1)

    for lst in delete_list:
        try:
            pfd_matrix.remove(lst)
        except ValueError:
            pass
        
    pfd_matrix = sorted(pfd_matrix, key=lambda pfd_matrix_entry: int(pfd_matrix_entry[3]))
    #print pfd_matrix
    domains = []
    for row in pfd_matrix:
        domains.append(row[-2]) #save the pfam domains for the .pfs file

    return pfd_matrix, domains                       
                        
@timeit
def write_pfs(pfs_handle, domains):
    for domain in domains:
        pfs_handle.write(domain+" ")
    pfs_handle.close()
    
@timeit
def write_pfd(pfd_handle, matrix):
    for row in matrix:
        row = "\t".join(row)
        pfd_handle.write(row+"\n")
        
    pfd_handle.close() 
    
@timeit
def writeout(handle, dct, keys):
    for key in keys:
        handle.write(key+"\n"+dct[key]+"\n")
    handle.close()

@timeit
def hmmscan(fastafile, outputdir, name, cores):

    #removed --noali par

    hmmscan_cmd = "hmmscan --cpu " + str(cores) + " --domtblout " + str(outputdir)\
     + "/" +str(name) + ".domtable --cut_tc Pfam-A.hmm " + str(fastafile) 
    if verbose == True:
        print hmmscan_cmd
        
    subprocess.check_output(hmmscan_cmd, shell=True)
    
@timeit
def get_domains(filename):
    handle = open(filename, 'r')
    domains = []
    for line in handle:
        if line[0] != "#":
            domains.append(filter(None, line.split(" "))[1])
            
    return domains

@timeit
def genbank_parser_hmmscan_call(gb_files, skip_hmmscan):
    """Extract the CDS from the antismash genbank clusters, and provide these coding regions to hmmscan"""
    
    gbk_group = {} #Will contain the gbk cluster as key, and the assigned group as a value together with the definition
    fasta_dict = {} #should make fasta headers more unique
    
    for gb_file in gb_files:
        
        CDS_keys = []
        outputbase = gb_file.split("/")[-1].replace(".gbk", "")
        outputfile = outputdir + "/" + outputbase + ".fasta"
        multifasta = open(outputfile, "w")

        gb_handle = open(gb_file, "r")
        
        #Parse the gbk file for the gbk_group dictionary
        in_cluster = False

        for line in gb_handle:
            if "DEFINITION  " in line:
                definition = line.strip().replace("DEFINITION  ", "")
            if "  cluster  " in line:
                in_cluster = True
                
            if "/product" in line and in_cluster == True:
                group= ""
                #group = line.strip().split(" ")[-1].replace("/product=", "").replace(" ", "").replace("\"", "")
                group = line.strip().split("=")[-1].replace("\"", "")
                print outputbase, group
                gbk_group[outputbase] = [group, definition]
                gb_handle.close()
                break
        
        try:
            gb_record = SeqIO.read(open(gb_file, "r"), "genbank")
        except ValueError:
            print "cannot find file", gb_file
        
            
        features = get_all_features_of_type(gb_record, "CDS")
        
        
        for feature in features:
            gene_id = ""
            protein_id = ""
            try:
                gene_id =  feature.qualifiers['gene']
            except KeyError:
                pass
            
            try:
                protein_id =  feature.qualifiers['protein_id']
            except KeyError:
                pass
            
            try:
                fasta_header = ">" + "loc:" + str(feature.location)\
                + ":gid:" + str(gene_id )\
                + ":pid:" + str(protein_id)\
                + ":loc_tag:" + str(feature.qualifiers['locus_tag'])
            except KeyError:
                print "no locus tag available in gbk file", gb_file
                fasta_header = ">" + "loc:" + str(feature.location)\
                + ":gid:" + str(gene_id )\
                + ":pid:" + str(protein_id)
                
            CDS_keys.append(fasta_header)
            fasta_dict[fasta_header] =\
            str(feature.qualifiers['translation'][0])
            

        writeout(multifasta, fasta_dict, CDS_keys)
        multifasta.close()
        
        if skip_hmmscan == False:
            hmmscan(outputfile, outputdir, outputbase, cores) #Run hmmscan
        else:
            print "Skipping hmmscan"
        
    return gbk_group, fasta_dict
            
@timeit       
def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

@timeit
def get_hmm_output_files():
    hmm_table_list = []
    for dirpath, dirnames, filenames in os.walk(str(outputdir) + "/"):
        for fname in filenames:

            if fname.split(".")[-1] == "domtable":
                if open(outputdir + "/" + fname, "r").readlines()[3][0] != "#": #if this is false, hmmscan has not found any domains in the sequence
                    hmm_table_list.append(fname)
                    if verbose == True:
                        print fname

                
                    
    return hmm_table_list


@timeit
def get_gbk_files(gbksamples):
    """Find .gbk files, and store the .gbk files in lists, separated by sample."""
    genbankfiles=[] #Will contain lists of gbk files
    dirpath = ""
    file_counter = 0
    for dirpath, dirnames, filenames in os.walk(str(os.getcwd()) + "/"):
        genbankfilelist=[]
        
        for fname in filenames:
            if fname.split(".")[-1] == "gbk" and "final" not in fname:
                file_counter += 1
                genbankfilelist.append(dirpath + "/" + fname)
                if verbose == True:
                    print fname
                if gbksamples == True:
                    genbankfiles.append(genbankfilelist)
                    genbankfilelist = []
                
        if genbankfilelist != []:
            genbankfiles.append(genbankfilelist)
    
    return genbankfiles

@timeit
def domtable_parser(gbk, hmm_table):
    """Parses the domain table output files from hmmscan"""
    
##example from domain table output:
# target name        accession   tlen query name                                    accession   qlen   E-value  score  bias   #  of  c-Evalue  i-Evalue  score  bias  from    to  from    to  from    to  acc description of target
#------------------- ---------- -----                          -------------------- ---------- ----- --------- ------ ----- --- --- --------- --------- ------ ----- ----- ----- ----- ----- ----- ----- ---- ---------------------
#Lycopene_cycl        PF05834.8    378 loc:[0:960](-):gid::pid::loc_tag:['ctg363_1'] -            320   3.1e-38  131.7   0.0   1   1   1.1e-40   1.8e-36  126.0   0.0     7   285    33   295    31   312 0.87 Lycopene cyclase protein


    pfd_matrix = []

    handle = open(hmm_table, 'r')
    for line in handle:
        if line[0] != "#":
            try:
                splitline = filter(None, line.split(" "))
                pfd_row = []
                pfd_row.append(gbk)         #add clustername or gbk filename
                
                pfd_row.append(splitline[13]) #add the score
    
    ##example of header_list ['loc', '[2341', '3538](+)', 'gid', '', 'pid', '', 'loc_tag', "['ctg4508_7"]]
    
                header_list = splitline[3].split(":")
                pfd_row.append(header_list[header_list.index("gid")+1]) #add gene ID if known
                pfd_row.append(splitline[19])#first coordinate, env coord from
                pfd_row.append(splitline[20])#second coordinate, env coord to
                loc_split = header_list[2].split("]") #second coordinate (of CDS) and the direction
                pfd_row.append(loc_split[1]) #add direction          
    
                pfd_row.append(splitline[1]) #pfam id
                pfd_row.append(splitline[0]) #domain name
                pfd_row.append(splitline[3])#cds header
                
                pfd_matrix.append(pfd_row)
            except ValueError:
                print "line: ", line
                print "file", hmm_table
            #print pfd_row

    return pfd_matrix

@timeit
def get_feature(string):
    feature = ""
    
    if re.search(r' {2,}\S* {2,}', string):
        feature = re.search(r' {2,}\S* {2,}', string).group(0)
        
    return feature

@timeit
def hmm_table_parser(gbk, hmm_table):
##AB050629.gbk	score   yxjC	1	1167	+	PF03600	CitMHS

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
##            pfd_row.append(splitline[-1])

##            row = "\t".join(pfd_row)
            
            pfd_matrix.append(pfd_row)

    return pfd_matrix

def get_anchor_domains(filename):
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
        

@timeit
def get_domain_list(filename):
    handle = open(filename, 'r')
    domains_string = handle.readline()
    domains = domains_string.split(" ")
    return domains

def make_domains_output_name(args):
    foldername = "domains_" + "_".join(args)
    return foldername.replace(" ", "_")  

def make_network_output_name(args):
    foldername = "networks_" + "_".join(args)
    return foldername.replace(" ", "_")


def CMD_parser():
    parser = OptionParser()
    
    
    parser.add_option("-o", "--outputdir", dest="outputdir", default="",
                      help="output directory, this contains your pfd,pfs,network and hmmscan output files")
    parser.add_option("-c", "--cores", dest="cores", default=8,
                      help="amount of cores hmmscan and the python script itself uses")
    parser.add_option("-l", "--limit", dest="limit", default=-1,
                      help="-limit- parameter of run_antismash")
    
    parser.add_option("-e", "--exclude_disc_nodes", dest="exclude_disc_nodes", action="store_false", default=True,
                      help="Exclude nodes that have no edges to other nodes from the network, default is true.")
    
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="toggle to true")
    parser.add_option("-d", "--domain_overlap_cutoff", dest="domain_overlap_cutoff", default=0.1,
                      help="Specify at which overlap percentage domains are considered to overlap")
    parser.add_option("--Jaccardw", dest="Jaccardw", default=0.4,
                      help="SJaccard weight")
    parser.add_option("--DDSw", dest="DDSw", default=0.5,
                      help="DDS weight")
    parser.add_option("--GKw", dest="GKw", default=0.1,
                      help="GK weight")
    parser.add_option("--domainsout", dest="domainsout", default="domains",
                      help="outputfolder of the pfam domain fasta files")
    parser.add_option("--mafft_pars", dest="mafft_pars", default="",
                      help="Add single/multiple parameters for mafft specific enclosed by quotation marks e.g. \"--nofft --parttree\"")
    parser.add_option("--al_method", dest="al_method", default="--retree 2",
                      help="alignment method for mafft, if there's a space in the method's name, enclose by quotation marks. default: \"--retree 1\"")
    parser.add_option("--maxiterate", dest="maxit", default=10,
                      help="maxiterate parameter in mafft, default is 0")
    parser.add_option("--mafft_threads", dest="mafft_threads", default=-1,
                      help="Set the number of threads in mafft, -1 sets the number of threads as the number of physical cores")

    parser.add_option("--use_perc_id", dest="use_perc_id", action="store_true", default=False,
                      help="Let the script calculate the percent identity between sequences? \
                      Or use the distout scores from the mafft output? - default")

    parser.add_option("--skip_hmmscan", dest="skip_hmmscan", action="store_true", default=False,
                      help="When skipping hmmscan, the GBK files should be available, and the domain tables need to be in the output folder.")

    parser.add_option("--sim_cutoffs", dest="sim_cutoffs", default="1,0.7,0.65,0.6,0.5,0.4,0.3,0.2,0.1",
                      help="Generate networks using multiple simmilarity (raw distance) cutoff values, example: \"1,0.5,0.1\"")
    
    parser.add_option("-a", "--anchorweight", dest="anchorweight", default=0.2,
                      help="Defines how much a pfam domain that's an anchor will contribute to the DDS of that distance.")
    
    parser.add_option("-n", "--nbhood", dest="nbhood", default=4,
                      help="")
    parser.add_option("-s", "--gbksamples", dest="gbksamples", action="store_true", default=False,
                      help="If each seperate gbk file represents a different sample, toggle to true.")
    

    (options, args) = parser.parse_args()
    return options, args



def gbk_processing(parms):
    gbks = parms[0]
    clusters = []
    #Loop over the samples
    #samplefolder = "/".join(gbks[0].split("/")[0:-1])
    samplename = ".".join(gbks[0].split("/")[-1].split(".")[0:-2])
    print "running hmmscan and or parsing the hmmscan output files on sample:", samplename 
    #outputdir = samplefolder + "/" + str(options.outputdir)         
    

    group_dct, fasta_dict = genbank_parser_hmmscan_call(gbks, skip_hmmscan) #runs hammscan and returns the CDS in the cluster
    hmm_domtables = get_hmm_output_files() 

    hmms = [] 
    for hmm_file in hmm_domtables:
        if samplename in hmm_file:
            hmms.append(hmm_file.replace(".domtable", "")) 
    
    clusters.append(hmms) #remember the clusters per sample        
    #loop over clusters
    for outputbase in hmms:
        pfsoutput = outputdir + "/" + outputbase + ".pfs"
        pfs_handle = open(pfsoutput, 'w')
        
        
        #pfd_matrix = hmm_table_parser(outputbase+".gbk", outputdir +"/"+ hmm_file)
        pfd_matrix = domtable_parser(outputbase, outputdir + "/" + outputbase+".domtable")
        filtered_matrix, domains = check_overlap(pfd_matrix, domain_overlap_cutoff)  #removes overlapping domains, and keeps the highest scoring domain
        save_domain_seqs(filtered_matrix, fasta_dict, domainsout) #save the sequences for the found domains per pfam domain
        write_pfs(pfs_handle, domains)
        
        BGC = BGC_dic_gen(filtered_matrix)
        BGCs[outputbase] = BGC
        pfdoutput = outputdir + "/" + outputbase + ".pfd"
        pfd_handle = open(pfdoutput, 'w')
        write_pfd(pfd_handle, filtered_matrix)
    

     
    #=======================================================================
    # network_matrix, networkfilename = generate_network(hmms, group_dct, str(samplename), networks_folder, "domain_dist", anchor_domains)
    # for cutoff in cutoff_list:
    #     write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", exclude_disc_nodes)
    #=======================================================================
    
    return group_dct.items(), clusters

@timeit
def main():
    options, args = CMD_parser()
    
    #anchor_handle = open("anchor_doms.txt", 'r') 
    print "Program starts"
    
    if options.outputdir == "":
        print "please provide a name for an output folder using parameter -o or --outputdir"
        sys.exit(0)
    
    anchor_domains = get_anchor_domains("anchor_domains.txt")
    
    global verbose
    global BGCs
    global DMS
    global timings_file
    global Jaccardw
    global DDSw
    global GKw
    global anchorweight
    global nbhood
    global cores
    global skip_hmmscan
    global domain_overlap_cutoff
    global domainsout
    global outputdir
    global group_dct
    exclude_disc_nodes = options.exclude_disc_nodes
    domain_overlap_cutoff = options.domain_overlap_cutoff
    skip_hmmscan = options.skip_hmmscan
    cores = int(options.cores)
    gbksamples = options.gbksamples
    nbhood = int(options.nbhood)
    anchorweight = float(options.anchorweight)
    Jaccardw = float(options.Jaccardw)
    DDSw = float(options.DDSw) 
    GKw = float(options.GKw)
    cutoff_list = options.sim_cutoffs.split(",")
    outputdir = str(options.outputdir)
    verbose = options.verbose
    args = sys.argv[1:]
    args.remove("-o")
    args.remove(outputdir)
    networks_folder = make_network_output_name(args)
    domainsout = make_domains_output_name(args)
    
    
    #print "$$$$$$$$$$$$$$", outputdir
    #subprocess.check_output("mkdir " + outputdir, shell=True)
    try:
        subprocess.check_output("mkdir " + outputdir, shell=True)
    except subprocess.CalledProcessError, e:
        pass
    
    #print "$$$$$$$$$$$$$$", outputdir + "/" + networks_folder
    #subprocess.check_output("mkdir " + outputdir + "/" + networks_folder, shell=True)
    try:
        subprocess.check_output("mkdir " + outputdir + "/" + networks_folder, shell=True)
    except subprocess.CalledProcessError, e:
        pass
    
    #print "$$$$$$$$$$$$$$", outputdir + "/" + domainsout
    #subprocess.check_output("mkdir " + outputdir + "/" + domainsout, shell=True)
    try:
        subprocess.check_output("mkdir " + outputdir + "/" + domainsout, shell=True)
    except subprocess.CalledProcessError, e:
        pass
    
    
    timings_file = open(outputdir + "/" + "runtimes.txt", 'w')
    print "getting gbk files"
    gbk_files = get_gbk_files(gbksamples) #files will contain lists of gbk files per sample. Thus a matrix contains lists with gbk files by sample.
    
    if gbk_files == []:
        print "No .gbk files were found"
    
    print gbk_files
    
    """BGCs -- dictionary of this structure:  BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""
     
     
    BGCs = {} #will contain the BGCs
    group_dct = {}
    clusters = [] #remember the clusters/bgcs per sample 
    
    parms=[]
    for gbks in gbk_files:
        parms.append([gbks])

    print "Running hmmscan and parsing the gbk files"
    pool = multiprocessing.Pool(cores) #create the appropriate amount of pool instances
    results = pool.map(gbk_processing, parms)   #Assigns the data to the different workers and pools the results
                                                #back into the network_matrix variable   
    
    
    for unpooled_result in results:
        group_dct = dict(group_dct.items() + unpooled_result[0]) #get the returned items, and merge them with the group_dct's items
        clusters += unpooled_result[1]
    
    print '$$$$$$$$$$$$$$$$$$$$', group_dct
    
    #parse group_dct and clusters from results

     
    print "Generating all_vs_all network with domain_dist"
    network_matrix, networkfilename = generate_network(list(iterFlatten(clusters)), group_dct, "all_vs_all", networks_folder, "domain_dist", anchor_domains)   
    for cutoff in cutoff_list:
        write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", exclude_disc_nodes)
        
        
        
    if verbose == True:
        print "BGCs:", BGCs
        BGC_handle = open("BGCs.txt", "w")
        for item in BGCs.items():
            BGC_handle.write(str(item)+"\n")
        BGC_handle.close()
        
        
        
    """DMS -- dictionary of this structure: DMS = {'general_domain_name_x': { ('specific_domain_name_1',
    'specific_domain_name_2'): (sequence_identity, alignment_length), ... }   }
        - general_domain_name_x: as above
        - ('specific_domain_name_1', 'specific_domain_name_2'): pair of specific domains, sorted alphabetically
        - (sequence_identity, alignment_length): sequence identity and alignment length of the domain pair"""



    DMS = {}
    fasta_domains = get_domain_fastas(domainsout)
    #Fill the DMS variable by using all 'domains.fasta'  files
    for domain_file in fasta_domains:
        print "Running MAFFT and parsing the distout file for domain:", domain_file
        domain=domain_file.replace(".fasta", "")
        run_mafft(options.al_method, options.maxit, options.mafft_threads, options.mafft_pars, domain)
        
        
        if options.use_perc_id == True:
            print "using percent identiy to calculate cluster diversity"
            fasta_handle = open(domain + ".algn", 'r')
            fasta_dict = fasta_parser(fasta_handle) #overwrites the fasta dictionary for each fasta file
           
            spec_domains_dict = {}
            for spec_domain in fasta_dict.keys():
                for spec_domain_nest in fasta_dict.keys():
                    if spec_domain != spec_domain_nest:
                        #tuple(sorted([seta[0],setb[0]]))
                        dist, length = calc_perc_identity(fasta_dict[spec_domain], fasta_dict[spec_domain_nest])
                        spec_domains_dict[tuple(sorted([spec_domain, spec_domain_nest]))] = (dist, length)
                    
            DMS[domain.split("/")[-1]] = spec_domains_dict
            
        else:      
            print "using distout file from mafft to calculate cluster diversity with sequence sim"
            domain_pairs = distout_parser(domain+".fasta" + ".hat2")
            if domain_pairs != {}:
                DMS[domain.split("/")[-1]] = domain_pairs
        
        
        
    #===========================================================================
    # #Compare the gene clusters within one sample, and save them in tab delimited .txt files.
    # for clusters_per_sample in clusters:
    #     network_matrix, networkfilename = generate_network(clusters_per_sample, group_dct, "", networks_folder, "seqdist", anchor_domains)
    #     for cutoff in cutoff_list:
    #         write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", exclude_disc_nodes)
    #===========================================================================
    

    
    network_matrix, networkfilename = generate_network(list(iterFlatten(clusters)), group_dct, "all_vs_all", networks_folder, "seqdist", anchor_domains)
    for cutoff in cutoff_list:
        write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", exclude_disc_nodes)
        
        
    if verbose == True:
        print "Saving the DMS variable to DMS.txt"
        DMS_handle = open("DMS.txt", "w")
        for item in DMS.items():
            DMS_handle.write(str(item)+"\n")
        DMS_handle.close()
    
  
main()
timings_file.close()
