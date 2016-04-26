#!/usr/bin/env python


"""
Student/programmer: Marley Yeong
marleyyeong@live.nl
supervisor: Marnix Medema
dependencies: hmmer, biopython, mafft, munkres, numpy

Usage:   place this script in the same folder as the antismash output files
         include the munkres.py file in this folder
         include the files necessary for hmmscan in this folder
         Example runstring: python ~/bgc_networks/bgc_networks.py -o exclude_small_bgcs --mafft_threads 12 -m 5000 -c 12 --include_disc_nodes

Status: development/testing

Todo: 
Calculate and report on some general properties of the generated networks

"""


from functions import *
import fileinput, pickle, sys, math
from munkres import Munkres
import numpy as np
import sys
import re
import os
import subprocess
from optparse import OptionParser
import math
import urllib2 
import multiprocessing
from functools import partial
import time
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
        insignificant_runtime = 1 #prevents an overload 
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        runtime = time2-time1
        
        runtime_string = '%s function took %0.3f s' % (f.func_name, runtime)
        
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
    
    cluster_file1 = output_folder + "/" + cluster1 + ".pfs"
    cluster_file2 = output_folder + "/" + cluster2 + ".pfs"
    
    
    A = get_domain_list(cluster_file1)
    B = get_domain_list(cluster_file2)
    
    
    if A == [''] or B == ['']:
        print "Regarding files", cluster_file1, cluster_file2
        print "One or both of these clusters contain no pfam domains"
        sys.exit()
        return [cluster1,cluster2,'','','','',str(1),str(1),str(1),'','']
    

    if dist_method == "seqdist":
        dist = cluster_distance(cluster1, cluster2, A, B, anchor_domains) #sequence dist
    elif dist_method == "domain_dist":
        dist = Distance_modified(A, B, 0, 4) #domain dist
        
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
def generate_network(bgc_list, group_dct, networkfilename, networks_folder, dist_method, anchor_domains, cores):
    #Contents of the network file: clustername1 clustername2 group1, group2, -log2score, dist, squared similarity
    "saves the distances as the log2 of the similarity"

    
    if networkfilename == "":
    
        networkfilename = output_folder + "/" + networks_folder + "/" + "networkfile_" + dist_method + "_" + \
        str("".join(bgc_list[0].split(".")[0:-2]))
    else:
        networkfilename = output_folder + "/" + networks_folder + "/" + "networkfile_" + dist_method + "_" + \
        networkfilename
        
    print "Generating network file with name:", networkfilename

    cluster_queue = []
    cluster_pairs = []
    
    for bgc in bgc_list:
        cluster_queue.append(bgc) #deep copy
        
    for cluster1 in bgc_list:
        cluster_queue = remove_values_from_list(cluster_queue, cluster1) # remove ALL instances of the bgc
        
        for cluster2 in cluster_queue:
            #addclusterpair
            cluster_pairs.append([cluster1, cluster2, dist_method, anchor_domains])    
            
    
    #maxtasksperchild is the number of tasks a worker process can complete before it will exit and be replaced with a fresh worker process, to
    #enable unused resources to be freed. The default maxtasksperchild is None, which means worker processes will live as long as the pool.
    
    pool = multiprocessing.Pool(cores, maxtasksperchild=1000) #create the appropriate amount of pool instances, with a limited amount of tasks per child process
    network_matrix = pool.map(generate_dist_matrix, cluster_pairs)  #Assigns the data to the different workers and pools the results
                                                                    #back into the network_matrix variable   

    sorted_network_matrix = sorted(network_matrix, key=lambda network_matrix_entry: network_matrix_entry[5]) #sort the matrix on the log2 scores
    return sorted_network_matrix, networkfilename
    
    
def cluster_distance(A, B, A_domlist, B_domlist, anchor_domains): 
  
    #key is name of GC, values is list of specific pfam domain names
    clusterA = BGCs[A] #will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    clusterB = BGCs[B]
     
    try:
        #calculates the intersect
        Jaccard = len(set(clusterA.keys()) & set(clusterB.keys())) / \
              float( len(set(clusterA.keys())) + len(set(clusterB.keys())) \
              - len(set(clusterA.keys()) & set(clusterB.keys())))
    except ZeroDivisionError:
        print "Zerodivisionerror during the Jaccard distance calculation. Can only happen when one or more clusters contains no pfam domains."
        print "keys of clusterA", A, clusterA.keys()
        print "keys of clusterB", B, clusterB.keys()
         
    intersect = set(clusterA.keys() ).intersection(clusterB.keys()) #shared pfam domains
    not_intersect = []
    
    #dom_seq_dist: Difference in sequence per domain. If one cluster doesn't have a domain at all, but the other does, 
    #this is a sequence difference of 1. If both clusters contain the domain once, and the sequence is the same, there is a seq diff of 0.
    #S: Max occurence of each domain
    dom_seq_dist_anchor,Sa = 0,0 #DDSa
    dom_seq_dist,S = 0,0 #dom_seq_dist
    pair = ""


    #start calculating the DDS


    for domain in set(clusterA.keys() + clusterB.keys()):
        if domain not in intersect:
            not_intersect.append(domain)
            
    for unshared_domain in not_intersect: #no need to look at seq identity or anchors, since these domains are unshared
        #for each occurence of an unshared domain do dom_seq_dist += count of domain and S += count of domain
        unshared_occurrences = []
        try:
            unshared_occurrences = clusterA[unshared_domain]
        except KeyError:
            unshared_occurrences = clusterB[unshared_domain]
            
        dom_seq_dist += len(unshared_occurrences)
        S += len(unshared_occurrences)
        
        
    for shared_domain in intersect:
        seta = clusterA[shared_domain]
        setb = clusterB[shared_domain]
        
        if len(seta+setb) == 2: #The domain occurs only once in both clusters   
            pair = tuple(sorted([seta[0],setb[0]]))
            
            try:
                seq_dist = 1-DMS[shared_domain][pair][0] #1-sequence_similarity
            except KeyError:
                print "KeyError on", pair
                
                errorhandle = open(str(A)+".txt", 'w')
                errorhandle.write(str(pair)+"\n")
                errorhandle.write(str(DMS[shared_domain]))
                errorhandle.close()
            
            if shared_domain.split(".")[0] in anchor_domains: 
                Sa += max(len(seta),len(setb))
                dom_seq_dist_anchor += seq_dist 
            else:
                S += max(len(seta),len(setb))
                dom_seq_dist += seq_dist
        else:                   #The domain occurs more than once in both clusters
            accumulated_distance = 0
            
            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            for domsa in range(len(seta)):
                for domsb in range(domsa, len(setb)):
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    
                    try:
                        Similarity = DMS[shared_domain][pair][0]
                    except KeyError:
                        print "KeyError on", pair
                        errorhandle = open(str(B)+".txt", 'w')
                        errorhandle.write(str(pair)+"\n")
                        errorhandle.write(str(DMS[shared_domain]))
                        errorhandle.close()
                    
                    seq_dist = 1-Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist
                
            #Only use the best scoring pairs
            Hungarian = Munkres()
            #print "DistanceMatrix", DistanceMatrix
            BestIndexes = Hungarian.compute(DistanceMatrix)
            #print "BestIndexes", BestIndexes
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            #print "accumulated_distance", accumulated_distance
                         
            sum_seq_dist = (abs(len(seta)-len(setb)) + accumulated_distance)  #essentially 1-sim
            
            if shared_domain.split(".")[0] in anchor_domains: 
                Sa += max(len(seta),len(setb))
                dom_seq_dist_anchor += sum_seq_dist 
            else:
                S += max(len(seta),len(setb))
                dom_seq_dist += sum_seq_dist 

 
    #  calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A_domlist]
    Ar.reverse()
    GK = max([calculate_GK(A_domlist, B_domlist, nbhood), calculate_GK(Ar, B_domlist, nbhood)])
    
    print "GKGKGKGGKGKGKGKGKGKGKGKGKGK"
    print GK
    print B_domlist
    print A_domlist
    

    
     
    if dom_seq_dist_anchor != 0:
        dom_seq_dist_combined = (anchorweight * (dom_seq_dist_anchor / float(Sa))) + ((1 - anchorweight) * (dom_seq_dist / float(S)))   #Recalculate dom_seq_dist by giving preference to anchor domains
        DDS = dom_seq_dist_combined / float(S)
    else:
        DDS = dom_seq_dist / float(S) 
        
        
    print "DDSDDSDDSDDSDDS"
    print DDS
    
    
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
                print "KeyError on", pair
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
                        print "KeyError on", pair
                      
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
    
    print mafft_cmd
    subprocess.check_output(mafft_cmd, shell=True)


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
def genbank_parser_hmmscan_call(gb_files, outputdir, cores, gbk_group, skip_hmmscan):
    """Extract the CDS from the antismash genbank clusters, and provide these coding regions to hmmscan"""
    
    #gbk_group = {} #Will contain the gbk cluster as key, and the assigned group as a value together with the definition
    fasta_dict = {} #should make fasta headers more unique
    
    for gb_file in gb_files:
        
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
        
        feature_counter = 0
        for feature in features:
            feature_counter += 1
            
            start = feature.location.start
            end = feature.location.end
            
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

            fasta_header = outputbase + "_ORF" + str(feature_counter)+ ":gid:" + str(gene_id) + ":pid:" + str(protein_id) + ":loc:" + str(start) + ":" + str(end)
            fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this
            
            #===================================================================
            # try:
            #     fasta_header = "loc:" + str(feature.location)\
            #     + ":gid:" + str(gene_id)\
            #     + ":pid:" + str(protein_id)\
            #     + ":loc_tag:" + str(feature.qualifiers['locus_tag'])
            # except KeyError:
            #     print "no locus tag available in gbk file", gb_file
            #     fasta_header = "loc:" + str(feature.location)\
            #     + ":gid:" + str(gene_id )\
            #     + ":pid:" + str(protein_id)
            #===================================================================
                
                
            
            fasta_header = ">"+(fasta_header.replace(" ", "")) #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header
                
            sequence = str(feature.qualifiers['translation'][0]) #in the case of the translation there should be one and only one entry (entry zero)
            if sequence != "":
                fasta_dict[fasta_header] = sequence
            
        dct_writeout(multifasta, fasta_dict) #save the coding sequences in a fasta format
        multifasta.close()
        
        if skip_hmmscan == False:
            hmmscan(outputfile, outputdir, outputbase, cores) #Run hmmscan
        else:
            print "Skipping hmmscan"
        
    return gbk_group, fasta_dict
            

def CMD_parser():
    parser = OptionParser()
    
    parser.add_option("-o", "--outputdir", dest="outputdir", default="",
                      help="Output directory, this will contain your pfd, pfs, network and hmmscan output files.")
    parser.add_option("-i", "--inputdir", dest="inputdir", default="",
                      help="Input directory of gbk files, if left empty, all gbk files in current and lower directories will be used.")
    parser.add_option("-c", "--cores", dest="cores", default=8,
                      help="Set the amount of cores the script and hmmscan may use")
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true", default=False,
                      help="Toggle to true, will print and save some variables.")
    parser.add_option("--include_disc_nodes", dest="include_disc_nodes", action="store_true", default=False,
                      help="Toggle to true. Include nodes that have no edges to other nodes from the network, default is false.")
    parser.add_option("-d", "--domain_overlap_cutoff", dest="domain_overlap_cutoff", default=0.1,
                      help="Specify at which overlap percentage domains are considered to overlap")
    parser.add_option("-m", "--min_bgc_size", dest="min_bgc_size", default=0,
                      help="Provide the minimum size of a bgc in base pairs, default is 0bp")
    parser.add_option("--seqdist_networks", dest="seqdist_networks", default="A",
                      help="Mode A generates the all vs all networks with sequence distance. Mode S compares clusters within a sample.\
                       Sample input: \"S,A\" generates the samplewise and all vs all. Default is \"A\"")
    parser.add_option("--domaindist_networks", dest="domaindist_networks", default="A",
                      help="Mode A generates the all vs all networks with domain distance. Mode S compares clusters within a sample.\
                       Sample input: \"A\" only generates the all vs all. Default is \"A\"")
    parser.add_option("--Jaccardw", dest="Jaccardw", default=0.1,
                      help="SJaccard weight")
    parser.add_option("--DDSw", dest="DDSw", default=0.7,
                      help="DDS weight")
    parser.add_option("--GKw", dest="GKw", default=0.2,
                      help="GK weight")
    parser.add_option("--domainsout", dest="domainsout", default="domains",
                      help="outputfolder of the pfam domain fasta files")
    parser.add_option("--anchorfile", dest="anchorfile", default="anchor_domains.txt",
                      help="Provide a custom name for the anchor domains file, default is anchor_domains.txt.")
    parser.add_option("--exclude_gbk_str", dest="exclude_gbk_str", default="final",
                      help="If this string occurs in the gbk filename, this will not be used for the analysis. Best to just leave out these samples to begin with.")
    parser.add_option("--mafft_pars", dest="mafft_pars", default="",
                      help="Add single/multiple parameters for mafft specific enclosed by quotation marks e.g. \"--nofft --parttree\"")
    parser.add_option("--al_method", dest="al_method", default="--retree 2",
                      help="alignment method for mafft, if there's a space in the method's name, enclose by quotation marks. default: \"--retree 1\"")
    parser.add_option("--maxiterate", dest="maxit", default=10,
                      help="Maxiterate parameter in mafft, default is 10")
    parser.add_option("--mafft_threads", dest="mafft_threads", default=-1,
                      help="Set the number of threads in mafft, -1 sets the number of threads as the number of physical cores")
    parser.add_option("--use_mafft_distout", dest="use_perc_id", action="store_false", default=True,
                      help="Let the script calculate the percent identity between sequences? \
                      Or use the distout scores from the mafft output? As default it calculates the percent identity from the MSAs.")
    parser.add_option("--skip_hmmscan", dest="skip_hmmscan", action="store_true", default=False,
                      help="When skipping hmmscan, the GBK files should be available, and the domain tables need to be in the output folder.")
    parser.add_option("--sim_cutoffs", dest="sim_cutoffs", default="1,0.7,0.65,0.6,0.5,0.4,0.3,0.2,0.1",
                      help="Generate networks using multiple similarity (raw distance) cutoff values, example: \"1,0.5,0.1\"")
    parser.add_option("-a", "--anchorweight", dest="anchorweight", default=0.1,
                      help="Weight of the anchor domains in the DDS distance metric. Default is set to 0.1.")
    parser.add_option("-n", "--nbhood", dest="nbhood", default=4,
                      help="nbhood variable for the GK distance metric, default is set to 4.")
    parser.add_option("-s", "--gbksamples", dest="gbksamples", action="store_true", default=False,
                      help="If each seperate gbk file represents a different sample, toggle to true. \
                      Saves you the effort of having to place each gbk file in a different folder.")

    (options, args) = parser.parse_args()
    return options, args


@timeit
def main():
    options, args = CMD_parser()
    
    if options.outputdir == "":
        print "please provide a name for an output folder using parameter -o or --outputdir"
        sys.exit(0)
    
    anchor_domains = get_anchor_domains(options.anchorfile)
    
    global verbose
    global BGCs
    global DMS
    global output_folder
    global timings_file
    global Jaccardw
    global DDSw
    global GKw
    global anchorweight
    global nbhood
    global cores
    include_disc_nodes = options.include_disc_nodes
    cores = int(options.cores)
    gbksamples = options.gbksamples
    nbhood = int(options.nbhood)
    anchorweight = float(options.anchorweight)
    Jaccardw = float(options.Jaccardw)
    DDSw = float(options.DDSw) 
    GKw = float(options.GKw)
    cutoff_list = options.sim_cutoffs.split(",")
    output_folder = str(options.outputdir)
    verbose = options.verbose
    args = sys.argv[1:]
    args.remove("-o")
    args.remove(output_folder)
    networks_folder = make_network_output_name(args)
    domainsout = make_domains_output_name(args)
    seqdist_networks = options.seqdist_networks.split(",")
    domaindist_networks = options.domaindist_networks.split(",")
    
    
    try:
        subprocess.check_output("mkdir " + output_folder, shell=True)
    except subprocess.CalledProcessError, e:
        pass

    try:
        subprocess.check_output("mkdir " + output_folder + "/" + networks_folder, shell=True)
    except subprocess.CalledProcessError, e:
        pass

    try:
        subprocess.check_output("mkdir " + output_folder + "/" + domainsout, shell=True)
    except subprocess.CalledProcessError, e:
        pass

    timings_file = open(output_folder + "/" + "runtimes.txt", 'wa') #open the file that will contain the timed functions
    
    
    gbk_files = get_gbk_files(options.inputdir, gbksamples, int(options.min_bgc_size), options.exclude_gbk_str) #files will contain lists of gbk files per sample. Thus a matrix contains lists with gbk files by sample.
    check_data_integrity(gbk_files)

    
    """BGCs -- dictionary of this structure:  BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""
     
     
    BGCs = {} #will contain the BGCs
    global group_dct
    group_dct = {}
    clusters = []
    
    #Loop over the samples 
    for gbks in gbk_files:
        #samplefolder = "/".join(gbks[0].split("/")[0:-1])
        samplename = ".".join(gbks[0].split("/")[-1].split(".")[0:-2]) #gbk files should look something like samplename.clusternumber.gbk
        print "running hmmscan and or parsing the hmmscan output files on sample:", samplename 
        #outputdir = samplefolder + "/" + str(options.outputdir)         
        
        group_dct, fasta_dict = genbank_parser_hmmscan_call(gbks, output_folder, cores, group_dct, options.skip_hmmscan) #runs hammscan and returns the CDS in the cluster
        hmm_domtables = get_hmm_output_files(output_folder) 
    
        hmms = [] 
        for hmm_file in hmm_domtables:
            if samplename in hmm_file:
                hmms.append(hmm_file.replace(".domtable", "")) 
        
        clusters.append(hmms) #remember the clusters per sample        
        #loop over clusters
        for outputbase in hmms:
            print "Processing domtable file:", outputbase
            pfsoutput = output_folder + "/" + outputbase + ".pfs"
            pfs_handle = open(pfsoutput, 'w')
            
            #pfd_matrix = hmm_table_parser(outputbase+".gbk", output_folder +"/"+ hmm_file)
            pfd_matrix = domtable_parser(outputbase, output_folder + "/" + outputbase + ".domtable")
            
            #===================================================================
            # pfdoutput = output_folder + "/" + outputbase + ".pfd2"
            # pfd_handle = open(pfdoutput, 'w')
            # write_pfd(pfd_handle, pfd_matrix)
            #===================================================================
            
            
            filtered_matrix, domains = check_overlap(pfd_matrix, options.domain_overlap_cutoff)  #removes overlapping domains, and keeps the highest scoring domain
            save_domain_seqs(filtered_matrix, fasta_dict, domainsout, output_folder, outputbase) #save the sequences for the found domains per pfam domain
            write_pfs(pfs_handle, domains)
            BGC = BGC_dic_gen(filtered_matrix)
            BGCs[outputbase] = BGC
            pfdoutput = output_folder + "/" + outputbase + ".pfd"
            pfd_handle = open(pfdoutput, 'w')
            write_pfd(pfd_handle, filtered_matrix)
        
        
        if "S" in domaindist_networks:
            print "Generating pairwise networks with domain_dist"
            network_matrix, networkfilename = generate_network(hmms, group_dct, str(samplename), networks_folder, "domain_dist", anchor_domains, cores)
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", include_disc_nodes)
     
     
    if "A" in domaindist_networks:
        print "Generating all_vs_all network with domain_dist"
        network_matrix, networkfilename = generate_network(list(iterFlatten(clusters)), group_dct, "all_vs_all", networks_folder, "domain_dist", anchor_domains, cores)   
        for cutoff in cutoff_list:
            write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", include_disc_nodes)
         
        
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
    fasta_domains = get_domain_fastas(domainsout, output_folder)
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
                        sim, length = calc_perc_identity(fasta_dict[spec_domain], fasta_dict[spec_domain_nest], spec_domain, spec_domain_nest, domain)
                        spec_domains_dict[tuple(sorted([spec_domain.replace(">",""), spec_domain_nest.replace(">","")]))] = (sim, length)
                    
            DMS[domain.split("/")[-1]] = spec_domains_dict
            
        else:      
            print "using distout file from mafft to calculate cluster diversity with sequence sim"
            domain_pairs = distout_parser(domain+".fasta" + ".hat2")
            if domain_pairs != {}:
                DMS[domain.split("/")[-1]] = domain_pairs
                
    if verbose == True:
        print "Saving the DMS variable to DMS.txt"
        DMS_handle = open("DMS.txt", "w")
        for item in DMS.items():
            DMS_handle.write(str(item)+"\n")
        DMS_handle.close()
    
    
    if "S" in seqdist_networks:     
        #Compare the gene clusters within one sample, and save them in tab delimited .txt files.
        for clusters_per_sample in clusters:
            network_matrix, networkfilename = generate_network(clusters_per_sample, group_dct, "", networks_folder, "seqdist", anchor_domains, cores)
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", include_disc_nodes)
                

    if "A" in seqdist_networks:
        network_matrix, networkfilename = generate_network(list(iterFlatten(clusters)), group_dct, "all_vs_all", networks_folder, "seqdist", anchor_domains, cores)
        for cutoff in cutoff_list:
            write_network_matrix(network_matrix, cutoff, networkfilename + "_c" + cutoff + ".network", include_disc_nodes)
         
        

    
    
main()
timings_file.close()
