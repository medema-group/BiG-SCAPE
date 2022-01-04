from multiprocessing import Pool
from functools import partial

from src.big_scape.clustering import generate_dist_matrix

# @timeit
def generate_network(cluster_pairs, cores, cluster_names, bgc_class_names, domain_list,
                     output_folder, gene_domain_count, corebiosynthetic_position,
                     bgc_gene_orientation, bgc_class_weight, anchor_domains, bgcs, mode, bgc_info,
                     aligned_domain_sequences, verbose, domains_folder):
    """Distributes the distance calculation part
    cluster_pairs is a list of triads (cluster1_index, cluster2_index, BGC class)
    """

    pool = Pool(cores, maxtasksperchild=100)

    #Assigns the data to the different workers and pools the results back into
    # the network_matrix variable
    # TODO: reduce argument count
    partial_func = partial(generate_dist_matrix, cluster_names=cluster_names,
                           bgc_class_names=bgc_class_names, domain_list=domain_list,
                           output_folder=output_folder, gene_domain_count=gene_domain_count,
                           corebiosynthetic_position=corebiosynthetic_position,
                           bgc_gene_orientation=bgc_gene_orientation,
                           bgc_class_weight=bgc_class_weight,
                           anchor_domains=anchor_domains, bgcs=bgcs, mode=mode, bgc_info=bgc_info,
                           aligned_domain_sequences=aligned_domain_sequences, verbose=verbose,
                           domains_folder=domains_folder)
    network_matrix = pool.map(partial_func, cluster_pairs)

    # --- Serialized version of distance calculation ---
    # For the time being, use this if you have memory issues
    #network_matrix = []
    #for pair in cluster_pairs:
      #network_matrix.append(generate_dist_matrix(pair))

    return network_matrix

def write_network_matrix(matrix, cutoffs_filenames, include_singletons, cluster_names, bgc_info):
    """
    An entry in the distance matrix is currently (all floats):
      0         1       2      3      4    5    6    7    8        9    10    11        12
    clus1Idx clus2Idx  rawD  sqrtSim  Jac  DSS  AI rDSSna  rDSSa   S    Sa lcsStartA lcsStartB
        13        14
    seedLength reverse

    The final row in the network file is currently:
      0      1      2     3      4   5   6     7       8    9   10    11       12
    clus1  clus2  rawD  sqrtSim  J  DSS  AI  rDSSna  rDSSa  S   Sa  combGrp  ShrdGrp
    """

    #Open file handles for each cutoff
    networkfiles = {}
    cutoffs, filenames = zip(*cutoffs_filenames)
    headers = ["Clustername 1",
               "Clustername 2",
               "Raw distance",
               "Squared similarity",
               "Jaccard index",
               "DSS index",
               "Adjacency index",
               "raw DSS non-anchor",
               "raw DSS anchor",
               "Non-anchor domains",
               "Anchor domains",
               "Combined group",
               "Shared group"]
    for cutoff, filename in cutoffs_filenames:
        networkfiles[cutoff] = open(filename, "w")
        networkfiles[cutoff].write("\t".join(headers))
        networkfiles[cutoff].write("\n")

    #Dictionaries to keep track of connected nodes, to know which are singletons
    cluster_set_all = {}
    cluster_set_connected = {}
    for cutoff in cutoffs:
        cluster_set_all[cutoff] = set()
        cluster_set_connected[cutoff] = set()

    for matrix_entry in matrix:
        gc1 = cluster_names[int(matrix_entry[0])]
        gc2 = cluster_names[int(matrix_entry[1])]
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
            row.append(" - ".join(sorted([clus1group, clus2group])))
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
            cluster_set_all[cutoff].add(gc1)
            cluster_set_all[cutoff].add(gc2)

            if row[2] < cutoff:
                cluster_set_connected[cutoff].add(gc1)
                cluster_set_connected[cutoff].add(gc2)

                networkfiles[cutoff].write("\t".join(map(str, row)) + "\n")


    #Add the nodes without any edges, give them an edge to themselves with a distance of 0
    if include_singletons:
        for cutoff in cutoffs:
            for gcs in cluster_set_all[cutoff] - cluster_set_connected[cutoff]:
                #Arbitrary numbers for S and Sa domains: 1 of each (logical would be 0,0 but
                # that could mess re-analysis with divisions-by-zero;
                numbers = [gcs, gcs, "0", "1", "1", "1", "1", "0", "0", "1", "1", "", ""]
                networkfiles[cutoff].write("\t".join(numbers) + "\n")

    #Close all files
    for networkfile in networkfiles.values():
        networkfile.close()
