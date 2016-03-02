#!/usr/bin/env python

"""
Student/programmer: Marley Yeong
marleyyeong@live.nl
supervisor: Marnix Medema

Usage: score a network

"""
from optparse import OptionParser
import sys
import os
import math


def get_network_files():
    network_files = []
    for dirpath, dirnames, filenames in os.walk(str(os.getcwd()) + "/"):
        for fname in filenames:

            if fname.split(".")[-1] == "network":
                try:
                    open(dirpath + "/" + fname, "r").readlines()[1]
                    network_files.append(dirpath + "/" + fname)
                except IndexError:
                    pass #means the network file is empty
                    
    return network_files

def normalize(OldMin, OldMax, NewMin, NewMax):
    
    OldRange = (max(values) - min(values))
    if OldRange == 0:
        NewValue = NewMin
        
    else:
        NewRange = (NewMax - NewMin)  
        NewValue = (((OldValue - OldMin) * NewRange) / OldRange) + NewMin
    
    return NewValue

def load_groups(filename):
    """Parse the network groups"""
    
    group_handle = open(filename, "r")
    groups_dict = {}
    number_of_groups = 0
    for line in group_handle:
        number_of_groups += 1
        group = line.split(":")[0]
        clusters = line.split(":")[1].strip().split(",")
        for cluster in clusters:
            groups_dict[cluster.lower()] = group

    return groups_dict, number_of_groups


def get_disconnected_subnetworks(filename):
    network_handle = open(filename, "r")
    network_handle.readline() #skip the header
    
    network_connections = {} #holds the disconnected subnetworks

    for line in network_handle:
        c1,c2 = line.split("\t")[0:2]
        try:
            if c2 not in network_connections[c1]:
                network_connections[c1].append(c2)
        except KeyError:
            network_connections[c1] = [c2]
        
        try:
            if c1 not in network_connections[c2]:
                network_connections[c2].append(c1)
        except KeyError:        
            network_connections[c2] = [c1]
    
    
    disconnected_subnetworks = []    
    
    genes = network_connections.keys()   
    for gene_cluster in genes:
        subnetwork = []
        subnetwork = [gene_cluster] + network_connections[gene_cluster]
        for node in subnetwork:
            if node != gene_cluster:
                
                for possible_outer_node in network_connections[node]:
                    if possible_outer_node not in subnetwork:
                        subnetwork.append(possible_outer_node) 
                del genes[genes.index(node)]    #no need to delete the node we are currently looping over in the first loop, 
                                                #its edges have already been implemented in the subnetwork
                                                #But we do need to delete nodes we are using to look up 'outer' edges
                                                #otherwise we'll loop over these nodes again, and they will make another subnetwork 
                                                #even though they have already been implemented in a subnetwork/cluster
                        
        disconnected_subnetworks.append(subnetwork)
     
    return len(disconnected_subnetworks)


def get_network_score(groups, number_of_groups, filename, clusters_threshold):
    #columns in network file: clustername1    clustername2    group1    definition    group2    definition    -log2score    raw distance    squared similarity    combined group    shared group
    network_handle = open(filename, "r")
    network_handle.readline() #skip the header
    
    within_distances =0
    between_distances =0
    within_edges = 0
    between_edges = 0
    cluster_occurence = []
    subnetworks = get_disconnected_subnetworks(filename) #returns the amount of disconnected subnetworks
    
    for line in network_handle:
        cluster_list = []
        group_list = []
        linelist = line.split("\t")
        cluster_list.append(linelist[3].lower().replace(" biosynthetic gene cluster",""))
        cluster_list.append(linelist[5].lower().replace(" biosynthetic gene cluster",""))
            
        for cluster in cluster_list:
            try:
                group = groups[cluster]
            except KeyError:
                for key in groups.keys():
                    if key in cluster:
                        #there is a space in the cluster definition
                        group = key
                        

            group_list.append(group)
            if cluster.lower() not in cluster_occurence and cluster != '':
                cluster_occurence.append(cluster.lower())
                        
        if group_list[0] == group_list[1]:
            within_distances += float(line[7])
            within_edges += 1
        else:
            between_distances += float(line[7])
            between_edges += 1
            
    network_handle.close()
    

    #===========================================================================
    # if between_edges > number_of_groups:
    #     between_edges -= number_of_groups # we do want some between edges, mainly between the clusters, so we need to correct for this
    #     #might want to remove this
    # else:
    #     between_edges = 1
    #===========================================================================
    
    if between_edges == 0:
        between_edges = 1
    
    score = within_edges/float(between_edges)       
                                                        
    #===========================================================================
    # 
    # cluster_fraction = float(len(cluster_occurence)) / len(groups.keys())
    # score = float(score) * (cluster_fraction**5)            #apply penalty for the amount of lacking clusters
    #===========================================================================
    
    if subnetworks > number_of_groups:
        score = score / (((subnetworks + 1) - number_of_groups)**3)   #If more subnetworks exist than the amount groups
                                                                #apply a penalty
    try:
        return filename.split("/")[-2] + "\t" + filename.split("/")[-1] + "\t" + str(score)  + "\t" + str(subnetworks)
    except IndexError:
        return filename + "\t" + filename + "\t" + str(score) +  "\t" + str(subnetworks)
    

def CMD_parser():
    parser = OptionParser()
    parser.add_option("-n", "--networkfile", dest="networkfile", default="",
                      help="name of networkfile")
 
 
    parser.add_option("-m", "--clusters_threshold", dest="clusters_threshold", default=0.6,
                      help="name of networkfile")
 
    (options, args) = parser.parse_args()
    return options, args
 
 
if __name__=="__main__":
    
    options, args = CMD_parser()
    networkfile = options.networkfile
    scores_handle = open("network_scores.txt","w")
    groups,number_of_groups = load_groups("groups.txt")
    if options.networkfile == "":
        for network_file in get_network_files():
            score = get_network_score(groups, number_of_groups, network_file, float(options.clusters_threshold))
            if score != 0:
                scores_handle.write(score+"\n")
    else:
        print get_network_score(groups, number_of_groups, options.networkfile, float(options.clusters_threshold))
                
