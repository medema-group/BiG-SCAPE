import os
import sys

from collections import defaultdict
from itertools import combinations
from itertools import product as combinations_product

import networkx as nx

from src.big_scape.clustering import clusterJsonBatch
from src.big_scape.distance import write_distance_matrix, get_cluster_cache_async, get_cluster_cache, gen_dist_matrix_async, gen_dist_matrix
from src.big_scape.bgc_info import BgcInfo
from src.bgctools import sort_bgc
from src.utility import create_directory

def get_output_cutoffs_filenames(run, path_base, bgc_class):
    file_names = []
    for cutoff in run.cluster.cutoff_list:
        file_names.append(os.path.join(path_base, "{}_c{:.2f}.network".format(bgc_class, cutoff)))
    cutoffs_filenames = list(zip(run.cluster.cutoff_list, file_names))
    del file_names[:]
    return cutoffs_filenames


def reduce_network(network_matrix):
    reduced_network = []
    pos_alignments = {}
    for row in network_matrix:
        reduced_network.append([int(row[0]), int(row[1]), row[2]])
        reverse = False
        if row[-1] == 1.0:
            reverse = True
        pos_alignment = pos_alignments.setdefault(int(row[0]), {})
        # lcsStartA, lcsStartB, seedLength, reverse={True,False}
        pos_alignment[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
    return reduced_network, pos_alignments

def create_working_set(run, cluster_names, domain_list, bgc_info: BgcInfo, mix) -> dict:
    bgc_classes = defaultdict(list)

    if mix:
        bgc_classes["mix"] = []

    for cluster_idx, cluster_name in enumerate(cluster_names):
        if run.has_includelist:
            # extra processing because pfs info includes model version
            bgc_domain_set = set({x.split(".")[0] for x in domain_list[cluster_name]})

            if len(run.domain_includelist & bgc_domain_set) == 0:
                continue
        product = bgc_info[cluster_name].product
        predicted_class = sort_bgc(product)

        if predicted_class.lower() in run.valid_classes:
            if mix:
                bgc_classes["mix"].append(cluster_idx)
            else:
                bgc_classes[predicted_class].append(cluster_idx)

        # possibly add hybrids to 'pure' classes
        if not mix and run.options.hybrids:
            if predicted_class == "PKS-NRP_Hybrids":
                if "nrps" in run.valid_classes:
                    bgc_classes["NRPS"].append(cluster_idx)
                if "t1pks" in product and "pksi" in run.valid_classes:
                    bgc_classes["PKSI"].append(cluster_idx)
                if "t1pks" not in product and "pksother" in run.valid_classes:
                    bgc_classes["PKSother"].append(cluster_idx)

            if predicted_class == "Others" and "." in product:
                subclasses = set()
                for subproduct in product.split("."):
                    subclass = sort_bgc(subproduct)
                    if subclass.lower() in run.valid_classes:
                        subclasses.add(subclass)

                # Prevent mixed BGCs with sub-Others annotations to get
                # added twice (e.g. indole-cf_fatty_acid has already gone
                # to Others at this point)
                if "Others" in subclasses:
                    subclasses.remove("Others")


                for subclass in subclasses:
                    bgc_classes[subclass].append(cluster_idx)
                subclasses.clear()

    return bgc_classes


def generate_network(run, cluster_names, domain_list, bgc_info: BgcInfo, query_bgc, domain_count_per_gene,
                     corebiosynthetic_pos, bgc_gene_orientation, bgcs, aligned_domain_seqs,
                     mibig_set_indices, mibig_set, rundata_networks_per_run,
                     html_subs_per_run, mix=False):
    print("\n Working for each BGC class")

    # we have to find the idx of query_bgc
    if run.directories.has_query_bgc:
        if query_bgc in cluster_names:
            query_bgc_idx = cluster_names.index(query_bgc)
        else:
            sys.exit("Error finding the index of Query BGC")

    # we will need this for the two distance matrix generation calls

    cluster_cache = get_cluster_cache_async(run, cluster_names, domain_list, domain_count_per_gene, corebiosynthetic_pos, bgc_gene_orientation, bgcs.bgc_dict, bgc_info)

    # create working set with indices of valid clusters
    bgc_classes = create_working_set(run, cluster_names, domain_list, bgc_info, mix)

    if mix:
        print("\n  {} ({} BGCs)".format("Mix", str(len(bgc_classes["mix"]))))

        # create output directory
        create_directory(os.path.join(run.directories.network, "mix"), "  Mix", False)
    else:
        # Preparing gene cluster classes
        print("  Sorting the input BGCs\n")

        class_names_len = len(run.distance.bgc_class_names)
        bgc_class_name_2_index = dict(zip(run.distance.bgc_class_names, range(class_names_len)))

    # only make folders for the bgc_classes that are found
    for bgc_class in bgc_classes:
        if run.directories.has_query_bgc:
            # not interested in this class if our Query BGC is not here...
            if query_bgc_idx not in bgc_classes[bgc_class]:
                continue

        print("\n  {} ({} BGCs)".format(bgc_class, str(len(bgc_classes[bgc_class]))))
        if run.mibig.use_mibig:
            if len(set(bgc_classes[bgc_class]) & mibig_set_indices) == len(bgc_classes[bgc_class]):
                print(" - All clusters in this class are MIBiG clusters -")
                print("  If you'd like to analyze MIBiG clusters, turn off the --mibig option")
                print("  and point --inputdir to the Annotated_MIBiG_reference folder")
                continue

        # create output directory
        create_directory(os.path.join(run.directories.network, bgc_class), "  All - " + bgc_class, False)

        # Create an additional file with the final list of all clusters in the class
        print("   Writing annotation files")
        network_annotation_path = os.path.join(run.directories.network, bgc_class, "Network_Annotations_" + bgc_class + ".tsv")
        with open(network_annotation_path, "w") as network_annotation_file:
            network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
            for idx in bgc_classes[bgc_class]:
                bgc = cluster_names[idx]
                product = bgc_info[bgc].product
                network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")

        print("   Calculating all pairwise distances")
        if run.directories.has_query_bgc:
            pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx], bgc_classes[bgc_class])])
        else:
        # convert into a set of ordered tuples
            pairs = set([tuple(sorted(combo)) for combo in combinations(bgc_classes[bgc_class], 2)])

        if mix:
            cluster_pairs = [(x, y, -1) for (x, y) in pairs]
        else:
            cluster_pairs = [(x, y, bgc_class_name_2_index[bgc_class]) for (x, y) in pairs]

        pairs.clear()

        network_matrix = gen_dist_matrix_async(run, cluster_pairs,
                                               cluster_names,
                                               cluster_cache,
                                               aligned_domain_seqs)
        #pickle.dump(network_matrix,open("others.ntwrk",'wb'))
        del cluster_pairs[:]
        #network_matrix = pickle.load(open("others.ntwrk", "rb"))

        # add the rest of the edges in the "Query network"
        if run.directories.has_query_bgc:
            new_set = []

            # rows from the distance matrix that will be pruned
            del_list = []

            for idx, row in enumerate(network_matrix):
                bgc_a, bgc_b, distance = int(row[0]), int(row[1]), row[2]

                # avoid QBGC-QBGC
                if bgc_a == bgc_b:
                    continue

                if distance <= run.cluster.max_cutoff:
                    if bgc_a == query_bgc_idx:
                        new_set.append(bgc_b)
                    else:
                        new_set.append(bgc_a)
                else:
                    del_list.append(idx)

            for idx in sorted(del_list, reverse=True):
                del network_matrix[idx]
            del del_list[:]

            pairs = set([tuple(sorted(combo)) for combo in combinations(new_set, 2)])

            if mix:
                cluster_pairs = [(x, y, -1) for (x, y) in pairs]
            else:
                cluster_pairs = [(x, y, bgc_class_name_2_index[bgc_class]) for (x, y) in pairs]

            pairs.clear()
            network_matrix_new_set = gen_dist_matrix_async(run, cluster_pairs, cluster_names,
                                                           cluster_cache, aligned_domain_seqs)
            del cluster_pairs[:]

            # Update the network matrix (QBGC-vs-all) with the distances of
            # QBGC's GCF
            network_matrix.extend(network_matrix_new_set)

            # Update actual list of BGCs that we'll use
            bgc_classes[bgc_class] = new_set
            bgc_classes[bgc_class].extend([query_bgc_idx])
            bgc_classes[bgc_class].sort()

            # Create an additional file with the list of all clusters in the class + other info
            # This version of the file only has information on the BGCs connected to Query BGC
            print("   Writing annotation file (Query BGC)")
            network_annotation_path = os.path.join(run.directories.network, bgc_class, "Network_Annotations_" + bgc_class + "_QueryBGC.tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in bgc_classes[bgc_class]:
                    bgc = cluster_names[idx]
                    product = bgc_info[bgc].product
                    network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
        elif run.mibig.use_mibig:
            nx_graph = nx.Graph()
            nx_graph.add_nodes_from(bgc_classes[bgc_class])
            mibig_set_del = []
            network_matrix_set_del = []

            for idx, row in enumerate(network_matrix):
                bgc_a, bgc_b, distance = int(row[0]), int(row[1]), row[2]
                if distance <= run.cluster.max_cutoff:
                    nx_graph.add_edge(bgc_a, bgc_b, index=idx)

            for component in nx.connected_components(nx_graph): # note: 'component' is a set
                num_bgcs_subgraph = len(component)

                # catch if the subnetwork is comprised only of MIBiG BGCs
                if len(component & mibig_set_indices) == num_bgcs_subgraph:
                    for bgc in component:
                        mibig_set_del.append(bgc)

            # Get all edges between bgcs marked for deletion
            for (bgc_a, bgc_b, idx) in nx_graph.subgraph(mibig_set_del).edges.data('index'):
                network_matrix_set_del.append(idx)

            # delete all edges between marked bgcs
            for row_idx in sorted(network_matrix_set_del, reverse=True):
                del network_matrix[row_idx]
            del network_matrix_set_del[:]

            print("   Removing {} non-relevant MIBiG BGCs".format(len(mibig_set_del)))
            bgc_to_class_idx = {}
            for idx, bgc in enumerate(bgc_classes[bgc_class]):
                bgc_to_class_idx[bgc] = idx
            for bgc_idx in sorted(mibig_set_del, reverse=True):
                del bgc_classes[bgc_class][bgc_to_class_idx[bgc_idx]]
            del mibig_set_del[:]

        if mix:
            pass
        else:
            if len(bgc_classes[bgc_class]) < 2:
                continue

        print("   Writing output files")
        path_base = os.path.join(run.directories.network, bgc_class)
        cutoffs_filenames = get_output_cutoffs_filenames(run, path_base, bgc_class)
        write_distance_matrix(network_matrix, cutoffs_filenames, run.options.include_singletons, cluster_names, bgc_info)

        print("  Calling Gene Cluster Families")
        reduced_network, pos_alignments = reduce_network(network_matrix)

        family_data = clusterJsonBatch(bgc_classes[bgc_class], path_base, bgc_class,
            reduced_network, pos_alignments, cluster_names, bgc_info,
            mibig_set, run.directories.pfd, run.directories.bgc_fasta, domain_list,
            bgcs.bgc_dict, aligned_domain_seqs, domain_count_per_gene, bgc_gene_orientation,
            cutoffs=run.cluster.cutoff_list, clusterClans=run.options.clans, clanCutoff=run.options.clan_cutoff,
            htmlFolder=run.directories.network_html)
        for network_html_folder_cutoff in family_data:
            rundata_networks_per_run[network_html_folder_cutoff].append(family_data[network_html_folder_cutoff])
            if mix:
                pass
            else:
                if len(family_data[network_html_folder_cutoff]["families"]) > 0:
                    class_info = {"name" : bgc_class, "css" : bgc_class, "label" : bgc_class}
                    html_subs_per_run[network_html_folder_cutoff].append(class_info)
        del bgc_classes[bgc_class][:]
        del reduced_network[:]
