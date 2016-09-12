#!/usr/bin/env python


"""
Date: 16-05-2016
Student/programmer: Marley Yeong
marleyyeong@live.nl
Supervisor: Marnix Medema
Dependencies: hmmer, biopython, mafft, munkres, numpy

Usage:   place this script in the same folder as the antismash output files
         include the munkres.py file in this folder
         include the files necessary for hmmscan in this folder
         Make sure to delete/empty the 'domains' output folder if you are rerunning using the same parameters!
         Example runstring: python ~/bgc_networks/bgc_networks.py -o exclude_small_bgcs --mafft_threads 12 -m 5000 -c 12 --include_disc_nodes

Status: development/testing

Todo: 
Calculate and report on some general properties of the generated networks
"""

from functions import *
from itertools import chain,repeat
from multiprocessing import Pool,cpu_count
from glob import glob

from munkres import Munkres
import sys
import os
import subprocess
from optparse import OptionParser
import math


import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from site import setBEGINLIBPATH
import cPickle as pickle # for storing and retrieving dictionaries


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
    
    cluster_file1 = os.path.join(output_folder, cluster1 + ".pfs")
    cluster_file2 = os.path.join(output_folder, cluster2 + ".pfs")
    
    
    A = get_domain_list(cluster_file1)
    B = get_domain_list(cluster_file2)
    
    
    # this really shouldn't happen if we've filtered the gene cluster list already
    if A == [''] or B == ['']:
        print("   Warning: Regarding distance between clusters " + cluster1 + " and " + cluster2 + ":")
        if A == [""] and B == [""]:
            print("   None have identified domains. Distance cannot be calculated")
        if A == [""]:            
            print("   Cluster " + cluster1 + " has no identified domains. Distance cannot be calculated")
        else:
            print("   Cluster " + cluster2 + " has no identified domains. Distance cannot be calculated")
        return [cluster1, cluster2, group_dct[cluster1][0], group_dct[cluster1][1], \
            group_dct[cluster2][0],group_dct[cluster2][1], '1.0', '1.0', '0.0','1.0','1.0', '1.0']
    

    if dist_method == "seqdist":
        dist, jaccard, dds, gk = cluster_distance(cluster1, cluster2, A, B, anchor_domains) #sequence dist
    elif dist_method == "domain_dist":
        dist, jaccard, dds, gk = Distance_modified(A, B, 0, 4) #domain dist
        
    if dist == 0:
        logscore = float("inf")
    else:
        logscore = 0
        try:
            logscore = -math.log(dist, 2) #Write exception, ValueError
        except ValueError:
            print "calculating the logscore with distance", dist, "failed in function generate_network, using distance method", dist_method, networkfilename
            
    #clustername1 clustername2 group1, group2, -log2score, dist, squared similarity, j, dds, gk
    network_row = [str(cluster1), str(cluster2), group_dct[cluster1][0],group_dct[cluster1][1], \
        group_dct[cluster2][0],group_dct[cluster2][1], str(logscore), str(dist), str((1-dist)**2), \
        jaccard, dds, gk]
    
    return network_row
    
        
@timeit
def generate_network(bgc_list, group_dct, networks_folder, dist_method, anchor_domains, cores):
    #Contents of the network file: clustername1 clustername2 group1, group2, -log2score, dist, squared similarity
    "saves the distances as the log2 of the similarity"
    
    # select every different pair. Sort by name to avoid double calculations
    cluster_pairs = []
    for i in range(len(bgc_list)-1):
        for j in range(i+1, len(bgc_list)):
            if bgc_list[i] < bgc_list[j]:
                cluster_pairs.append([bgc_list[i], bgc_list[j], dist_method, anchor_domains])
            else:
                cluster_pairs.append([bgc_list[j], bgc_list[i], dist_method, anchor_domains])
    
    # --- Use multiprocessing to distribute distance calculation ---
    #maxtasksperchild is the number of tasks a worker process can complete 
    # before it will exit and be replaced with a fresh worker process, to
    # enable unused resources to be freed. The default maxtasksperchild is 
    # None, which means worker processes will live as long as the pool.
    
    #create the appropriate amount of pool instances, with a limited amount of
    # tasks per child process
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
        
    #sorted_network_matrix = sorted(network_matrix, key=lambda network_matrix_entry: network_matrix_entry[5]) #sort the matrix on the log2 scores
    
    return network_matrix_dict
    
    
       
def cluster_distance(A, B, A_domlist, B_domlist, anchor_domains): 
    """Compare two clusters using information on their domains, and the sequences of the domains"""    

  
    #key is name of GC, values is list of specific pfam domain names
    clusterA = BGCs[A] #will contain a dictionary where keys are pfam domains, and values are domains that map to a specific sequence in the DMS variable
    clusterB = BGCs[B]
    
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
    pair = "" #pair of clusters to access their sequence identity
    
    # Case 1
    for unshared_domain in not_intersect: #no need to look at seq identity, since these domains are unshared
        #for each occurence of an unshared domain do domain_difference += count of domain and S += count of domain
        unshared_occurrences = []
        try:
            unshared_occurrences = clusterA[unshared_domain]
        except KeyError:
            unshared_occurrences = clusterB[unshared_domain]
            
        # don't look at domain version, hence the split
        if unshared_domain.split(".")[0] in anchor_domains:
            domain_difference_anchor += len(unshared_occurrences)
        else:
            domain_difference += len(unshared_occurrences)
        
    S = domain_difference # can be done because it's the first use of these
    S_anchor = domain_difference_anchor
        
        
    for shared_domain in intersect:
        seta = clusterA[shared_domain]
        setb = clusterB[shared_domain]
        
        # Case 2: The shared domains occurs only once in each gene cluster
        if len(seta+setb) == 2: #The domain occurs only once in both clusters   
            pair = tuple(sorted([seta[0],setb[0]]))
            
            try:
                seq_dist = 1-DMS[shared_domain][pair][0] #1-sequence_similarity
            except KeyError:
                print("For some reason, a sequence similarity value could not be retrieved from DMS...")
                print(shared_domain + ": " + pair)
                print(str(DMS[shared_domain]))
                sys.exit("(Case 2)")
            
            if shared_domain.split(".")[0] in anchor_domains: 
                S_anchor += 1
                domain_difference_anchor += seq_dist 
            else:
                S += 1
                domain_difference += seq_dist
                
        # Case 3: The domain occurs more than once in both clusters
        else:
            accumulated_distance = 0
            
            # Fill distance matrix between domain's A and B versions
            DistanceMatrix = [[1 for col in range(len(setb))] for row in range(len(seta))]
            for domsa in range(len(seta)):
                for domsb in range(len(setb)):
                    pair = tuple(sorted([seta[domsa], setb[domsb]]))
                    
                    try:
                        Similarity = DMS[shared_domain][pair][0] # the [1] key holds the sequence length
                    except KeyError:
                        print("For some reason, a sequence similarity value could not be retrieved from DMS...")
                        print(shared_domain + ": " + pair)
                        print(str(DMS[shared_domain]))
                        sys.exit("(Case 3)")
                        
                    seq_dist = 1-Similarity
                    DistanceMatrix[domsa][domsb] = seq_dist
                
            #Only use the best scoring pairs
            Hungarian = Munkres()
            #print "DistanceMatrix", DistanceMatrix
            BestIndexes = Hungarian.compute(DistanceMatrix)
            #print "BestIndexes", BestIndexes
            accumulated_distance = sum([DistanceMatrix[bi[0]][bi[1]] for bi in BestIndexes])
            #print "accumulated_distance", accumulated_distance
            
            # the difference in number of domains accounts for the "lost" (or not duplicated) domains
            sum_seq_dist = (abs(len(seta)-len(setb)) + accumulated_distance)  #essentially 1-sim
            
            if shared_domain.split(".")[0] in anchor_domains: 
                S_anchor += max(len(seta),len(setb))
                domain_difference_anchor += sum_seq_dist
            else:
                S += max(len(seta),len(setb))
                domain_difference += sum_seq_dist

        
    if domain_difference_anchor != 0 and domain_difference != 0:
        DDS = (anchorweight * (domain_difference_anchor / float(S_anchor))) + ((1 - anchorweight) * (domain_difference / float(S)))   #Recalculate domain_difference by giving preference to anchor domains
    elif domain_difference_anchor == 0:
        DDS = domain_difference / float(S)
    else: #only anchor domains were found
        DDS = domain_difference_anchor / float(S_anchor)
 
 
    # GK INDEX
    #  calculate the Goodman-Kruskal gamma index
    Ar = [item for item in A_domlist]
    Ar.reverse()
    GK = max([calculate_GK(A_domlist, B_domlist, nbhood), calculate_GK(Ar, B_domlist, nbhood)])
        
    
    DDS = 1-DDS #transform into similarity
    Distance = 1 - (Jaccardw * Jaccard) - (DDSw * DDS) - (GKw * GK) 
    if Distance < 0:
        print("Negative distance detected!")
        print("J: " + str(Jaccard) + ", DDS: " + str(DDS) + ", GK: " + str(GK))
        print("Jw: " + str(Jaccardw) + ", DDSw: " + str(DDSw) + ", GKw: " + str(GKw))
        sys.exit()
        
    return Distance, Jaccard, DDS, GK



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
        print(mafft_cmd)
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
        print "Generating fasta for: ", outputbase
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

def runHmmScan(fastaPath,hmmPath,outputdir):
    ## will run hmmscan command on a fasta file with a single core and generate a domtable file
    hmmFile = os.path.join(hmmPath,"Pfam-A.hmm")
    if os.path.isfile(fastaPath):
        name = fastaPath.split(os.sep)[-1].replace(".fasta","")
        outputName = os.path.join(outputdir, name+".domtable")
        if os.path.isfile(hmmFile):
            hmmscan_cmd = "hmmscan --cpu 1 --domtblout %s --cut_tc %s %s" % (outputName,hmmFile,fastaPath)
            if verbose == True:
                print("\t"+hmmscan_cmd)
            subprocess.check_output(hmmscan_cmd, shell=True)
        else:
            "Error: invalid path to HMM"
    else:
        "Error: Fasta File Doesn't Exist"

def parseHmmScan(hmmscanResults,outputdir,overlapCutoff):
    outputbase = hmmscanResults.split(os.sep)[-1].replace(".domtable", "")
    # try to read the domtable file to find out if this gbk has domains. Needed to parse domains into fastas anyway.
    if os.path.isfile(hmmscanResults):
        pfd_matrix = domtable_parser(outputbase, hmmscanResults)
        num_domains = len(pfd_matrix) # these might still be overlapped, but we only need a minimum number

        if num_domains > 0:
            print(" Processing domtable file: " + outputbase)

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
        print "Error hmmscan file doesn't exist"

def genbank_grp_dict(gb_files,gbk_group):
    for gb_file in gb_files:
        outputbase = gb_file.split(os.sep)[-1].replace(".gbk", "")
        gb_handle = open(gb_file, "r")

        #Parse the gbk file for the gbk_group dictionary
        # AntiSMASH-produced genbank files use the "cluster" tag
        # to annotate the Gene Cluster class in the "product" subtag...
        # (There should be just one "cluster" tag per genbank file,
        # if not, maybe you're using the wrong file with the whole seq.)
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
                #print(" " + outputbase + " " + group)
                gbk_group[outputbase] = [group, definition]
                gb_handle.close()
                break

        #...but other genbank files (e.g. MiBIG) do NOT have this tag
        if not in_cluster:
            gbk_group[outputbase] = ["no type", definition]
            gb_handle.close()
    return gbk_group

@timeit
def genbank_parser_hmmscan_call(gb_files, outputdir, cores, gbk_group, skip_hmmscan):
    """Extract the CDS from the antismash genbank clusters, and provide these coding regions to hmmscan"""
    
    #gbk_group = {} #Will contain the gbk cluster as key, and the assigned group as a value together with the definition
    fasta_dict = {} #should make fasta headers more unique
    
    for gb_file in gb_files:        
        sample_dict = {}
        
        outputbase = gb_file.split(os.sep)[-1].replace(".gbk", "")
        outputfile = os.path.join(outputdir, outputbase + ".fasta")
        gb_handle = open(gb_file, "r")
        
        #Parse the gbk file for the gbk_group dictionary
        # AntiSMASH-produced genbank files use the "cluster" tag
        # to annotate the Gene Cluster class in the "product" subtag...
        # (There should be just one "cluster" tag per genbank file, 
        # if not, maybe you're using the wrong file with the whole seq.)
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
                #print(" " + outputbase + " " + group)
                gbk_group[outputbase] = [group, definition]
                gb_handle.close()
                break
            
        #...but other genbank files (e.g. MiBIG) do NOT have this tag
        if not in_cluster:
            gbk_group[outputbase] = ["no type", definition]
            gb_handle.close()
        
        # write fasta sequence of features and predict domains
        if not skip_hmmscan:
            multifasta = open(outputfile, "w")
            
            # possible errors were catched previously in check_data_integrity
            features = get_all_features_of_type(gb_file, "CDS")
            
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
                
                strand = "+" if feature.strand == 1 else "-"

                fasta_header = outputbase + "_ORF" + str(feature_counter)+ ":gid:" + str(gene_id) + ":pid:" + str(protein_id) + ":loc:" + str(start) + ":" + str(end) + ":strand:" + strand
                fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this

                fasta_header = ">"+(fasta_header.replace(" ", "")) #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header
                    
                sequence = str(feature.qualifiers['translation'][0]) #in the case of the translation there should be one and only one entry (entry zero)
                if sequence != "":
                    fasta_dict[fasta_header] = sequence
                    sample_dict[fasta_header] = sequence #this should be used for hmmscan
                    

                
            dct_writeout(multifasta, sample_dict) #save the coding sequences in a fasta format
            multifasta.close()
        
        
            hmmscan(pfam_dir, outputfile, outputdir, outputbase, cores) #Run hmmscan
        
    return gbk_group, fasta_dict
            

def CMD_parser():
    parser = OptionParser()
    
    parser.add_option("-o", "--outputdir", dest="outputdir", default="",
                      help="Output directory, this will contain your pfd, pfs, network and hmmscan output files.")
    parser.add_option("-i", "--inputdir", dest="inputdir", default="",
                      help="Input directory of gbk files, if left empty, all gbk files in current and lower directories will be used.")
    parser.add_option("-c", "--cores", dest="cores", default=cpu_count(),
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
    
    parser.add_option("--domaindist_networks", dest="domaindist_networks", default="",
                      help="Mode A generates the all vs all networks with domain distance. Mode S compares clusters within a sample.\
                       Sample input: \"A\" only generates the all vs all. Default is \"\"")
    
    parser.add_option("--Jaccardw", dest="Jaccardw", default=0.2,
                      help="SJaccard weight, default is 0.2")
    parser.add_option("--DDSw", dest="DDSw", default=0.75,
                      help="DDS weight, default is 0.75")
    parser.add_option("--GKw", dest="GKw", default=0.05,
                      help="GK weight, default is 0.05")
    parser.add_option("-a", "--anchorw", dest="anchorweight", default=0.1,
                      help="Weight of the anchor domains in the DDS distance metric. Default is set to 0.1.")
    
    #parser.add_option("--domainsout", dest="domainsout", default="domains",
                      #help="outputfolder of the pfam domain fasta files")
    parser.add_option("--pfam_dir", dest="pfam_dir",
                      default=os.path.dirname(os.path.realpath(__file__)), 
                      help="Location of hmmpress-processed Pfam files. Default is same location of BiG-SCAPE")
    parser.add_option("--anchorfile", dest="anchorfile", default="anchor_domains.txt",
                      help="Provide a custom name for the anchor domains file, default is anchor_domains.txt.")
    parser.add_option("--exclude_gbk_str", dest="exclude_gbk_str", default="final",
                      help="If this string occurs in the gbk filename, this will not be used for the analysis. Best to just leave out these samples to begin with.")
    
    parser.add_option("--mafft_pars", dest="mafft_pars", default="",
                      help="Add single/multiple parameters for mafft specific enclosed by quotation marks e.g. \"--nofft --parttree\"")
    parser.add_option("--al_method", dest="al_method", default="--retree 2",
                      help="alignment method for mafft, if there's a space in the method's name, enclose by quotation marks. default: \"--retree 2\" corresponds to the FFT-NS-2 method")
    parser.add_option("--maxiterate", dest="maxit", default=1000,
                      help="Maxiterate parameter in mafft, default is 1000, corresponds to the FFT-NS-2 method")
    parser.add_option("--mafft_threads", dest="mafft_threads", default=-1,
                      help="Set the number of threads in mafft, -1 sets the number of threads as the number of physical cores")
    parser.add_option("--use_mafft_distout", dest="use_perc_id", action="store_false", default=True,
                      help="Let the script calculate the percent identity between sequences? \
                      Or use the distout scores from the mafft output? As default it calculates the percent identity from the MSAs.")
    
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

    parser.add_option("--splitProcess",dest="splitProcess",action="store_true",default=False,
                      help="Experimental Method: Splits initial work into separate processes equal to the number of cores used. "
                           "Namely (1) parsing genbanks, (2) running hmmscan and generating the domtable files, "
                           "(3) Parsing the domtable files to generate the pfs and pfd files. Will use these files to generate the BGC and DMS dictionary to rejoin the main workflow.")
    (options, args) = parser.parse_args()
    return options, args


if __name__=="__main__":
    options, args = CMD_parser()
    
    if options.outputdir == "":
        print "please provide a name for an output folder using parameter -o or --outputdir"
        sys.exit(0)
    
    anchor_domains = get_anchor_domains(options.anchorfile)
    
    global verbose
    global BGCs
    global DMS
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
    if anchorweight > 1.0 or anchorweight < 0.0:
        sys.exit("Invalid anchorweight parameter (must be between 0 and 1)")
    Jaccardw = float(options.Jaccardw)
    DDSw = float(options.DDSw) 
    GKw = float(options.GKw)
    cutoff_list = options.sim_cutoffs.split(",")
    if "1" not in cutoff_list:
        cutoff_list.append("1") # compulsory for re-runs
        print("Adding cutoff=1.0 case by default")
    output_folder = str(options.outputdir)
    pfam_dir = str(options.pfam_dir)
    verbose = options.verbose
    args = sys.argv[1:]
    if "-o" in args:
        args.remove("-o")
    else:
        args.remove("--outputdir")
    args.remove(output_folder)
    networks_folder = "networks"
    domainsout = "domains"
    seqdist_networks = options.seqdist_networks.split(",")
    domaindist_networks = options.domaindist_networks.split(",")
    # variable "samples" takes care of structure of files in get_gbk_files
    samples = True if "S" in seqdist_networks or "S" in domaindist_networks else False
    
    if options.skip_all:
        if options.skip_hmmscan or options.skip_mafft:
            print("Overriding --skip_hmmscan/--skip_mafft with --skip_all parameter")
            options.skip_hmmscan = False
            options.skip_mafft = False
    
    
    # Obtain gbk files
    global gbk_files
    # gbk_files will contain lists of gbk files per sample. Thus a matrix contains lists with gbk files by sample.
    gbk_files, sample_name = get_gbk_files(options.inputdir, samples, int(options.min_bgc_size), options.exclude_gbk_str) 
    check_data_integrity(gbk_files)
    
    
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
                print(" Emptying domains directory first")
                for thing in os.listdir(os.path.join(output_folder, domainsout)):
                    os.remove(os.path.join(output_folder, domainsout, thing))
            else:
                print(" Using existing domains directory")
        else:
            print("Fatal error when trying to create domains' directory")
            sys.exit(str(e))
    print("")


    timings_file = open(os.path.join(output_folder, "runtimes.txt"), 'w') #open the file that will contain the timed functions
    
    
    #===========================================================================
    # gbk_handle = open("gbk.txt", "w")
    # for i in gbk_files:
    #     gbk_handle.write(str(i)+"\n")
    # gbk_handle.close()
    #===========================================================================

    
    """BGCs -- dictionary of this structure:  BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""
     
     
    BGCs = {} #will contain the BGCs
    global group_dct
    group_dct = {}
    clusters = [] # main structure with samples/gene clusters to be analized later on
    sequences_per_domain = {} # to avoid calling MAFFT if only 1 seq. in a particular domain
    
    #Loop over the samples
    if options.skip_hmmscan or options.skip_mafft or options.skip_all:
        print("Skipping processing of hmmscan output files")
    else:
        print("Running hmmscan and parsing the hmmscan output files...")


    if options.splitProcess:
        if verbose:
            print("Trying Threading on %i cores" % cores)
        # since the structure of gbk_files is nested lists of file paths need to unpack this

        # This will go into clusters later for checking if we dropped files in processing
        genbankFileNames = set(chain(*gbk_files))
        baseNames = set(x.split(os.sep)[-1].replace(".gbk","") for x in genbankFileNames)

        ### Step 1: Generate Fasta Files
        print "Parsing genbank Files to generate fasta files for hmmsearch"

        # filter through task list to avoid unecessary computation, if output file is there and nonempty exclude it from list

        alreadyDone = set()

        for genbank in genbankFileNames:
            outputbase  = genbank.split(os.sep)[-1].replace(".gbk","")
            outputfile = os.path.join(output_folder,outputbase + '.fasta')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(genbank)

        if len(alreadyDone) > 0:
            print "Found output files for %s, skipping these" % list(x.split(os.sep)[-1].replace(".gbk","") for x in alreadyDone)


        # Generate Pool of workers
        pool = Pool(cores,maxtasksperchild=32)
        for genbankFile in set(genbankFileNames) - alreadyDone:
            pool.apply_async(generateFasta,args =(genbankFile,output_folder))
        pool.close()
        pool.join()

        print "Done generating fasta files checking output...."

        ### Step 2: Run hmmscan
        fastaFiles = glob(os.path.join(output_folder,"*.fasta"))
        fastaBases = set(x.split(os.sep)[-1].replace(".fasta","") for x in fastaFiles)

        if len(baseNames - fastaBases) > 0:
            print "WARNING no fasta files for:", baseNames - fastaBases

        if fastaFiles:
            alreadyDone = set()
            for fasta in fastaFiles:
                outputbase  = fasta.split(os.sep)[-1].replace(".fasta","")
                outputfile = os.path.join(output_folder,outputbase + '.domtable')
                if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                    alreadyDone.add(fasta)

            if len(alreadyDone) > 0:
                print "Found output files for %s, skipping these" % list(x.split(os.sep)[-1].replace(".fasta","") for x in alreadyDone)

            pool = Pool(cores,maxtasksperchild=1)
            for fastaFile in set(fastaFiles) - alreadyDone:
                pool.apply_async(runHmmScan,args=(fastaFile,pfam_dir,output_folder))
            pool.close()
            pool.join()

            print "Done generating domtable files checking output...."
        else:
            "ERROR: Tasklist empty for running hmmscan"

        ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files

        domtblFiles = glob(os.path.join(output_folder,"*.domtable"))
        domtblBases = set(x.split(os.sep)[-1].replace(".domtable","") for x in domtblFiles)
        if len(baseNames - domtblBases) > 0:
            print "WARNING no domtbl files for:", baseNames - domtblBases

        if domtblFiles:

            alreadyDone = set()

            for domtable in domtblFiles:
                outputbase  = fasta.split(os.sep)[-1].replace(".domtable","")
                outputfile = os.path.join(output_folder,outputbase + '.pfd')
                if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                    alreadyDone.add(domtable)

            if len(alreadyDone) > 0:
                print "Found output files for %s, skipping these" % list(x.split(os.sep)[-1].replace(".fasta","") for x in alreadyDone)

            pool = Pool(cores,maxtasksperchild=32)

            for domtableFile in set(domtblFiles) - alreadyDone:
                pool.apply_async(parseHmmScan, args=(domtableFile,output_folder,options.domain_overlap_cutoff))
            pool.close()
            pool.join()

            print "Done generating generating pfs and pfd files checking output...."
        else:
            "ERROR: Tasklist empty for parsing hmmscan"

        ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per sample objects
            ## Note: this code currently doesn't modify the fasta files to include the domain fasta IDs

        pfdFiles = glob(os.path.join(output_folder,"*.pfd"))
        pfdBases = set(x.split(os.sep)[-1].replace(".pfd","") for x in pfdFiles)

        if len(baseNames - pfdBases) > 0:
            print "WARNING no domtbl files for:", baseNames - pfdBases

        for outputbase in pfdBases:
            print "Processing: %s" % outputbase

            pfdFile = os.path.join(output_folder, outputbase + ".pfd")
            pfsFile = os.path.join(output_folder, outputbase + ".pfs")
            fasta_file = os.path.join(output_folder, outputbase + ".fasta")

            fasta_dict = fasta_parser(open(fasta_file))
            filtered_matrix = [map(lambda x: x.strip(), line.split('\t')) for line in open(pfdFile)]

            save_domain_seqs(filtered_matrix, fasta_dict, domainsout, output_folder, outputbase)

            BGCs[outputbase] = BGC_dic_gen(filtered_matrix)


        sample_name_new = []
        for gbks,sampleName in zip(gbk_files,chain(sample_name,repeat([]))):
            baseNames = set(gbk.split(os.sep)[-1].replace(".gbk", "") for gbk in gbks)
            clusters_per_sample = list(baseNames & pfdBases)
            if clusters_per_sample:
                group_dct = genbank_grp_dict(gbks,group_dct)
                clusters.append(clusters_per_sample)
                sample_name_new.append(sampleName)
        sample_name = sample_name_new

    ### This is the current workflow to generate the BGC dictionary, pfs, and pfd files #################

    else:
        for gbks in gbk_files:
            # the "group" dictionary must always be obtained as it's used when
            # writing the network file, so this call has to stay despite
            # --skip_hmmscan, --skip_mafft or --skip_all
            group_dct, fasta_dict = genbank_parser_hmmscan_call(gbks, output_folder, cores, group_dct, (options.skip_hmmscan or options.skip_mafft or options.skip_all) ) #runs hammscan and returns the CDS in the cluster

            clusters_per_sample = []
            for gbk in gbks:
                outputbase = gbk.split(os.sep)[-1].replace(".gbk", "")

                # try to read the domtable file to find out if this gbk has domains. Needed to parse domains into fastas anyway.
                pfd_matrix = domtable_parser(outputbase, os.path.join(output_folder, outputbase + ".domtable"))
                num_domains = len(pfd_matrix) # these might still be overlapped, but we only need a minimum number

                if num_domains > 0:
                    clusters_per_sample.append(outputbase)
                else:
                    print(" Could not find domains in " + outputbase + ".domtable file. Removing it from further analysis")

                if num_domains > 0 and not (options.skip_hmmscan or options.skip_mafft or options.skip_all):
                    print(" Processing domtable file: " + outputbase)

                    #pfd_matrix = hmm_table_parser(outputbase+".gbk", output_folder +"/"+ hmm_file)
                    filtered_matrix, domains = check_overlap(pfd_matrix, options.domain_overlap_cutoff)  #removes overlapping domains, and keeps the highest scoring domain
                    save_domain_seqs(filtered_matrix, fasta_dict, domainsout, output_folder, outputbase) #save the sequences for the found domains per pfam domain

                    # Save list of domains per BGC
                    pfsoutput = os.path.join(output_folder, outputbase + ".pfs")
                    pfs_handle = open(pfsoutput, 'w')
                    write_pfs(pfs_handle, domains)

                    # Save more complete information of each domain per BGC
                    pfdoutput = os.path.join(output_folder, outputbase + ".pfd")
                    pfd_handle = open(pfdoutput, 'w')
                    write_pfd(pfd_handle, filtered_matrix)

                    BGCs[outputbase] = BGC_dic_gen(filtered_matrix)

                    # keep track of number of sequences per domain. If only 1
                    #  sequence, skip MAFFT call for this particular domain.
                    for row in filtered_matrix:
                        try:
                            sequences_per_domain[row[5]] += 1
                        except KeyError:
                            sequences_per_domain[row[5]] = 1

            # empty entire samples could arise from individual gbks with no domains
            if len(clusters_per_sample) > 0:
                clusters.append(clusters_per_sample)

        print("")
    ##################################################################################################
        
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
    
    
    # Dictionary with pairwise distance information
    network_matrix = {}
    network_matrix_sample = {}
    
    # Distance without taking sequence similarity between specific domains into account
    if "S" in domaindist_networks or "A" in domaindist_networks:
        if options.skip_all: #read already calculated distances
            network_matrix = network_parser(os.path.join(output_folder, networks_folder, "networkfile_domain_dist_all_vs_all_c1.network"), Jaccardw, DDSw, GKw)
        
        # If user wants all-vs-all, no need to recalculate anything for sample networks, so, this goes first
        if "A" in domaindist_networks:
            print("* Generating all-vs-all network with domain distance method")
            if not options.skip_all:
                print(" Calculating all pairwise distances")
                network_matrix = generate_network(list(iterFlatten(clusters)), group_dct, networks_folder, "domain_dist", anchor_domains, cores)   
            
            # write network files
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, os.path.join(output_folder, networks_folder, "networkfile_domain_dist_all_vs_all_c" + cutoff + ".network"), include_disc_nodes)
    

        if "S" in domaindist_networks:
            if len(clusters) == 1 and "A" in domaindist_networks:
                print("* NOT generating networks per sample (only one sample, covered in the all-vs-all case)")
            else:
                print("* Generating sample networks with domain distance method")
                sn = 0
                for clusters_per_sample in clusters:
                    print(" Sample: " + sample_name[sn])
                    if len(clusters_per_sample) == 1:
                        try:
                            print(" Warning: Sample size = 1 detected. Not generating networks for this sample (" + sample_name[sn] + ")")
                        except IndexError:
                            print(sn)
                            print(sample_name)
                            sys.exit()
                    else:
                        network_matrix_sample.clear()
                        if "A" in domaindist_networks or options.skip_all:
                            # we should have the distances already calculated either from the all-vs-all case or retrieved from the file
                            for i in range(len(clusters_per_sample)-1):
                                for j in range(i+1, len(clusters_per_sample)):
                                    network_matrix_sample[tuple(sorted( [clusters_per_sample[i], clusters_per_sample[j]] ))] = \
                                        network_matrix[tuple(sorted( [clusters_per_sample[i], clusters_per_sample[j]] ))]
                        else:
                            print("  Calculating distances for this sample")
                            network_matrix_sample = generate_network(clusters_per_sample, group_dct, networks_folder, "domain_dist", anchor_domains, cores)
                    
                        # write network files
                        for cutoff in cutoff_list:
                            write_network_matrix(network_matrix_sample, cutoff, os.path.join(output_folder, networks_folder, "networkfile_domain_dist_" + sample_name[sn] + "_c" + cutoff + ".network"), include_disc_nodes)
                    sn += 1

    
    # Check whether user wants seqdist method networks before calculating DMS
    if "A" in seqdist_networks or "S" in seqdist_networks:
        print("")
        if options.skip_all:
            network_matrix = network_parser(os.path.join(output_folder, networks_folder, "networkfile_seqdist_all_vs_all_c1.network"), Jaccardw, DDSw, GKw)
            
        elif options.skip_mafft:
            DMS = {}
            print("Trying to read domain alignments (DMS.dict file)")
            with open(os.path.join(output_folder, "DMS.dict"), "r") as DMS_file:
                DMS = pickle.load(DMS_file)
                DMS_file.close()
                
        else:
            DMS = {}
            fasta_domains = get_domain_fastas(domainsout, output_folder)
            if options.skip_hmmscan or options.splitProcess: # in this case we didn't have a chance to fill sequences_per_domain

                for domain_fasta in fasta_domains:
                    domain_name = domain_fasta.split(os.sep)[-1].replace(".fasta", "")
                    
                    # fill fasta_dict
                    fasta_handle = open(domain_fasta, "r")
                    fasta_dict = fasta_parser(fasta_handle)
                    fasta_handle.close()
                    
                    sequences_per_domain[domain_name] = len(fasta_dict)
                        
                        
            print("Aligning domains")
            if options.use_perc_id == True:
                print("(Using calculated percentage identity for cluster diversity)")
            else:
                print("(Using percentage identity from MAFFT for cluster diversity)")
                
                
            """DMS -- dictionary of this structure: DMS = {'general_domain_name_x': { ('specific_domain_name_1',
            'specific_domain_name_2'): (sequence_identity, alignment_length), ... }   }
                - general_domain_name_x: as above
                - ('specific_domain_name_1', 'specific_domain_name_2'): pair of specific domains, sorted alphabetically
                - (sequence_identity, alignment_length): sequence identity and alignment length of the domain pair"""
            
            #Fill the DMS variable by using all 'domains.fasta'  files
            for domain_file in fasta_domains:
                domain_name = domain_file.split(os.sep)[-1].replace(".fasta", "")
                domain = domain_file.replace(".fasta", "")
                
                # avoid calling MAFFT if it's not possible to align (only one sequence)
                if sequences_per_domain[domain_name] == 1:
                    if verbose:
                        print(" Skipping MAFFT for domain " + domain_name + " (only one sequence)")
                else:                
                    if verbose:
                        print(" Running MAFFT for domain: " + domain_name)
                    
                    run_mafft(options.al_method, options.maxit, options.mafft_threads, options.mafft_pars, domain)
                    
                    if options.use_perc_id == True:
                        fasta_handle = open(domain + ".algn", 'r')
                        fasta_dict = fasta_parser(fasta_handle) #overwrites the fasta dictionary for each fasta file
                        fasta_handle.close()
                    
                        spec_domains_dict = {}
                        keys = fasta_dict.keys()
                        for i in range(len(fasta_dict)-1):
                            for j in range(i+1, len(fasta_dict)):
                                spec_domain = keys[i]
                                spec_domain_nest = keys[j]

                                sim, length = calc_perc_identity(fasta_dict[spec_domain], fasta_dict[spec_domain_nest], spec_domain, spec_domain_nest, domain)
                                spec_domains_dict[tuple(sorted([spec_domain.replace(">",""), spec_domain_nest.replace(">","")]))] = (sim, length)
                                
                        DMS[domain.split(os.sep)[-1]] = spec_domains_dict
                        
                    else:      
                        domain_pairs = distout_parser(domain + ".fasta" + ".hat2")
                        if domain_pairs != {}:
                            DMS[domain.split(os.sep)[-1]] = domain_pairs
            
            # though we don't really get to load this file by the time being
            with open(os.path.join(output_folder, "DMS.dict"), "w") as DMS_file:
                pickle.dump(DMS, DMS_file)
                DMS_file.close()
            print("")
            
            
        if "A" in seqdist_networks:
            print("* Generating all-vs-all network with domain-sequence distance method")
            if not options.skip_all:
                print(" Calculating all pairwise distances")
                network_matrix = generate_network(list(iterFlatten(clusters)), group_dct, networks_folder, "seqdist", anchor_domains, cores)
            for cutoff in cutoff_list:
                write_network_matrix(network_matrix, cutoff, os.path.join(output_folder, networks_folder, "networkfile_seqdist_all_vs_all_c" + cutoff + ".network"), include_disc_nodes)
                
        if "S" in seqdist_networks:
            if len(clusters) == 1 and "A" in seqdist_networks:
                print("* NOT generating networks per sample (only one sample, covered in the all-vs-all case)")
            else:
                print("* Generating sample networks with domain-sequence distance method")
                sn = 0
                for clusters_per_sample in clusters:
                    print(" Sample: " + sample_name[sn])
                    if len(clusters_per_sample) == 1:
                        print(" Warning: Sample size = 1 detected. Not generating network for this sample (" + sample_name[sn] + ")")
                    else:
                        network_matrix_sample.clear()
                        
                        if "A" in seqdist_networks or options.skip_all:
                            for i in range(len(clusters_per_sample)-1):
                                for j in range(i+1, len(clusters_per_sample)):
                                    network_matrix_sample[tuple(sorted([clusters_per_sample[i], clusters_per_sample[j]]))] = \
                                        network_matrix[tuple(sorted([clusters_per_sample[i], clusters_per_sample[j]]))]
                        else:
                            print("  Calculating distances for this sample")
                            network_matrix_sample = generate_network(clusters_per_sample, group_dct, networks_folder, "seqdist", anchor_domains, cores)
                    
                            
                        for cutoff in cutoff_list:
                            write_network_matrix(network_matrix_sample, cutoff, os.path.join(output_folder, networks_folder, "networkfile_seqdist_" + sample_name[sn] + "_c" + cutoff + ".network"), include_disc_nodes)
                    sn += 1
    timings_file.close()