from multiprocessing import Pool
from functools import partial

from src.utility.profiling import timeit
from src.big_scape.clustering import generate_dist_matrix

# @timeit
def generate_network(cluster_pairs, cores, clusterNames, bgcClassNames, DomainList, output_folder, DomainCountGene, corebiosynthetic_position, BGCGeneOrientation,
    bgc_class_weight, anchor_domains, BGCs, mode, bgc_info,
    AlignedDomainSequences, verbose, domains_folder):
    """Distributes the distance calculation part
    cluster_pairs is a list of triads (cluster1_index, cluster2_index, BGC class)
    """
    
    pool = Pool(cores, maxtasksperchild=100)
    
    #Assigns the data to the different workers and pools the results back into
    # the network_matrix variable
    partial_func = partial(generate_dist_matrix, clusterNames=clusterNames, bgcClassNames=bgcClassNames,
    DomainList=DomainList, output_folder=output_folder, DomainCountGene=DomainCountGene,
    corebiosynthetic_position=corebiosynthetic_position, BGCGeneOrientation=BGCGeneOrientation,
    bgc_class_weight=bgc_class_weight, anchor_domains=anchor_domains, BGCs=BGCs, mode=mode, bgc_info=bgc_info,
    AlignedDomainSequences=AlignedDomainSequences, verbose=verbose, domains_folder=domains_folder)
    network_matrix = pool.map(partial_func, cluster_pairs)

    # --- Serialized version of distance calculation ---
    # For the time being, use this if you have memory issues
    #network_matrix = []
    #for pair in cluster_pairs:
      #network_matrix.append(generate_dist_matrix(pair))

    return network_matrix

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


