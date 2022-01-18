import os
from array import array

from multiprocessing import Pool
from functools import partial

from src.big_scape.cluster_info import ClusterInfo
from src.pfam.misc import get_domain_list
from src.big_scape.scores import calc_distance_lcs


def get_cluster_cache_async(run, cluster_names, domain_list, gene_domain_count,
corebiosynthetic_positions, bgc_gene_orientation, bgcs):
    pool = Pool(run.options.cores, maxtasksperchild=100)

    # there are a couple of things we can calculate in advance and just once
    # these will be contained in these cluster_info classes
    func_preprocess = partial(preprocess_cluster_info, run=run, domain_list=domain_list, gene_domain_count=gene_domain_count,
    corebiosynthetic_positions=corebiosynthetic_positions, bgc_gene_orientation=bgc_gene_orientation,
    bgcs=bgcs)
    cluster_cache_list: list(ClusterInfo) = pool.map(func_preprocess, cluster_names)
    cluster_cache = dict()
    for cluster_info in cluster_cache_list:
        cluster_cache[cluster_info.cluster_name] = cluster_info
    
    return cluster_cache

def preprocess_cluster_info(cluster_name, run, domain_list, gene_domain_count,
corebiosynthetic_positions, bgc_gene_orientation, bgcs):
    try:
        domlist = domain_list[cluster_name]
    except KeyError:
        print(" Warning: domain list for {} was not found. Extracting from pfs files".format(cluster_name))

        cluster_file1 = os.path.join(run.directories.output, cluster_name + ".pfs")

        domlist = get_domain_list(cluster_file1)

    # "Domain Count per Gene". List of simple labels (integers) indicating number
    # of domains belonging to each gene
    dcg = gene_domain_count[cluster_name]

    # Position of the anchor genes (i.e. genes with domains in the anchor
    # domain list). Should probably be the Core Biosynthetic genes marked by
    # antiSMASH
    core_pos = corebiosynthetic_positions[cluster_name]

    # go = "gene orientation"
    go = bgc_gene_orientation[cluster_name]

    cluster_info = ClusterInfo(cluster_name, domlist, dcg, core_pos, go)
    
    cluster_info.init_dom_seq_slices(bgcs)

    #initialize domlist borders for AI
    cluster_info.init_dom_borders()

    cluster_info.init_gene_string()

    return cluster_info

# @timeit
def gen_dist_matrix_async(run, cluster_pairs, cluster_names,
                          cluster_info_cache, gene_domain_count,
                          bgc_gene_orientation, bgcs, bgc_info,
                          aligned_domain_sequences):
    """Distributes the distance calculation part
    cluster_pairs is a list of triads (cluster1_index, cluster2_index, BGC class)
    """

    pool = Pool(run.options.cores, maxtasksperchild=100)

    #Assigns the data to the different workers and pools the results back into
    # the network_matrix variable
    # TODO: reduce argument count
    func_dist_matrix = partial(generate_dist_matrix, run=run, cluster_names=cluster_names,
                               cluster_info_cache=cluster_info_cache,
                               gene_domain_count=gene_domain_count,
                               bgc_gene_orientation=bgc_gene_orientation,
                               bgcs=bgcs, bgc_info=bgc_info,
                               aligned_domain_sequences=aligned_domain_sequences)
    network_matrix = pool.map(func_dist_matrix, cluster_pairs)

    # --- Serialized version of distance calculation ---
    # For the time being, use this if you have memory issues
    # network_matrix = []
    # for pair in cluster_pairs:
    #     network_matrix.append(generate_dist_matrix(pair, run, cluster_names=cluster_names,
    #                        domain_list=domain_list,
    #                        gene_domain_count=gene_domain_count,
    #                        corebiosynthetic_position=corebiosynthetic_position,
    #                        bgc_gene_orientation=bgc_gene_orientation,
    #                        bgcs=bgcs, bgc_info=bgc_info,
    #                        aligned_domain_sequences=aligned_domain_sequences))

    return network_matrix

def generate_dist_matrix(parms, run, cluster_names, cluster_info_cache,
                         gene_domain_count, bgc_gene_orientation,
                         bgcs, bgc_info, aligned_domain_sequences):
    """Unpack data to actually launch cluster_distance for one pair of BGCs"""

    cluster_1_idx, cluster_2_idx, bgc_class_idx = [int(parm) for parm in parms]

    bgc_class = run.distance.bgc_class_names[bgc_class_idx]

    cluster_name_a = cluster_names[cluster_1_idx]
    cluster_name_b = cluster_names[cluster_2_idx]

    cluster_info_a: ClusterInfo = cluster_info_cache[cluster_name_a]
    cluster_info_b: ClusterInfo = cluster_info_cache[cluster_name_b]

    # this really shouldn't happen if we've filtered domain-less gene clusters already
    if len(cluster_info_a.domlist) == 0 or len(cluster_info_b.domlist) == 0:
        print("   Warning: Regarding distance between clusters {} and {}:".format(cluster_name_a, cluster_name_b))
        if len(cluster_info_a.domlist) == 0 and len(cluster_info_b.domlist) == 0:
            print("   None have identified domains. Distance cannot be calculated")
        elif (cluster_info_a.domlist) == 0:
            print("   Cluster {} has no identified domains. Distance set to 1".format(cluster_name_a))
        else:
            print("   Cluster {} has no identified domains. Distance set to 1".format(cluster_name_b))

        # last two values (S, Sa) should really be zero but this could give rise to errors when parsing
        # the network file (unless we catched the case S = Sa = 0

        # cluster1Idx, cluster2Idx, distance, jaccard, DSS, AI, rDSSNa, rDSSa,
        #   S, Sa, lcsStartA, lcsStartB
        return array('f', [cluster_1_idx, cluster_2_idx, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0])

    weights = run.distance.bgc_class_weight[bgc_class]

    dist, jaccard, dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB, seedLength, reverse = calc_distance_lcs(run,
                cluster_info_a, cluster_info_b, weights, bgcs, gene_domain_count, bgc_gene_orientation,
                bgc_info, aligned_domain_sequences)

    network_row = array('f', [cluster_1_idx, cluster_2_idx, dist, (1-dist)**2, jaccard,
                              dss, ai, rDSSna, rDSS, S, Sa, lcsStartA, lcsStartB,
                              seedLength, reverse])
    return network_row



def write_distance_matrix(distance_matrix, cutoffs_filenames, include_singletons, cluster_names, bgc_info):
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
    distance_file = {}
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
        distance_file[cutoff] = open(filename, "w")
        distance_file[cutoff].write("\t".join(headers))
        distance_file[cutoff].write("\n")

    #Dictionaries to keep track of connected nodes, to know which are singletons
    cluster_set_all = {}
    cluster_set_connected = {}
    for cutoff in cutoffs:
        cluster_set_all[cutoff] = set()
        cluster_set_connected[cutoff] = set()

    for matrix_entry in distance_matrix:
        gene_cluster_a = cluster_names[int(matrix_entry[0])]
        gene_cluster_b = cluster_names[int(matrix_entry[1])]
        row = [gene_cluster_a, gene_cluster_b]

        # get AntiSMASH annotations
        cluster_group_a = bgc_info[gene_cluster_a].product
        cluster_group_b = bgc_info[gene_cluster_b].product

        # add all the other floats
        row.extend(matrix_entry[2:-6])

        # add number of anchor/non-anchor domains as integers
        row.append(int(matrix_entry[-6]))
        row.append(int(matrix_entry[-5]))

        # prepare combined group
        if cluster_group_a != "" and cluster_group_b != "": #group1, group2
            row.append(" - ".join(sorted([cluster_group_a, cluster_group_b])))
        elif cluster_group_b != "":
            row.append(cluster_group_b)
        elif cluster_group_a != "":
            row.append(cluster_group_a)
        else:
            row.append("NA")

        # prepare share group (if they indeed share it)
        if cluster_group_a == cluster_group_b:
            row.append(cluster_group_a)
        else:
            row.append("")

        for cutoff in cutoffs:
            cluster_set_all[cutoff].add(gene_cluster_a)
            cluster_set_all[cutoff].add(gene_cluster_b)

            if row[2] < cutoff:
                cluster_set_connected[cutoff].add(gene_cluster_a)
                cluster_set_connected[cutoff].add(gene_cluster_b)

                distance_file[cutoff].write("\t".join(map(str, row)) + "\n")


    #Add the nodes without any edges, give them an edge to themselves with a distance of 0
    if include_singletons:
        for cutoff in cutoffs:
            for gcs in cluster_set_all[cutoff] - cluster_set_connected[cutoff]:
                #Arbitrary numbers for S and Sa domains: 1 of each (logical would be 0,0 but
                # that could mess re-analysis with divisions-by-zero;
                numbers = [gcs, gcs, "0", "1", "1", "1", "1", "0", "0", "1", "1", "", ""]
                distance_file[cutoff].write("\t".join(numbers) + "\n")

    #Close all files
    for distance_file in distance_file.values():
        distance_file.close()
