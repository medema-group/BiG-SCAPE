import logging
import os
import sys

from collections import defaultdict
from itertools import combinations
from itertools import product as combinations_product
from array import array

import networkx as nx

import pandas as pd
import numpy as np

from src.big_scape.scores import calc_jaccard
from src.big_scape.bgc_collection import BgcCollection
from src.big_scape.clustering import cluster_json_batch
from src.big_scape.distance import write_distance_matrix, gen_dist_matrix_async
from src.legacy.bgctools import sort_bgc
from src.utility import create_directory
from src.data.functions import get_bgc_ids, get_bgc_names, get_hmm_ids, get_features
from src.big_scape.cosine import get_corr_cosine_dists

def get_output_cutoffs_filenames(run, path_base, bgc_class):
    """Generate filenames for cutoffs in the run details

    Inputs:
        run: run details for this execution of BiG-SCAPE
        path_base: base path for network files
        bgc_class: the class to generate filenames for
    """
    file_names = []
    for cutoff in run.cluster.cutoff_list:
        file_names.append(os.path.join(path_base, "{}_c{:.2f}.network".format(bgc_class, cutoff)))
    cutoffs_filenames = list(zip(run.cluster.cutoff_list, file_names))
    del file_names[:]
    return cutoffs_filenames


def reduce_network(network_matrix):
    """reduce network matrix to only the relevant fieldsfor clustering

    inputs: network_matrix. a list of lists in the shape of [pair [dist, ]]
    """
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

def create_working_set(run, bgc_collection: BgcCollection, mix) -> dict:
    """Creates a working set of bgc class dictionaries that will later be looped through in
    pairwise distance calculation

    Inputs:
        run: run details for this execution of BiG-SCAPE
        bgc_collection: BgcCollection object which contains the collection of BGCs which will
            be used in pairwise comparison
        mix: boolean indicating if an additional 'mix' class should be created which will contain
            all BGCs. Note: this means that the comparison space is essentially doubled, and the
            comparisons increase exponentially."""
    bgc_classes = defaultdict(list)

    if mix:
        bgc_classes["mix"] = []

    for cluster_idx, cluster_name in enumerate(bgc_collection.bgc_name_tuple):
        if run.has_includelist:
            # extra processing because pfs info includes model version
            bgc_domain_set = set({x.split(".")[0] for x in bgc_collection.bgc_ordered_domain_list[cluster_name]})

            if len(run.domain_includelist & bgc_domain_set) == 0:
                continue
        product = bgc_collection.bgc_collection_dict[cluster_name].bgc_info.product

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

def generate_unrelated_row(cluster_1_idx, cluster_2_idx,):
    return array('f', [cluster_1_idx, cluster_2_idx, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0])

def generate_network(run, database, bgc_collection: BgcCollection, aligned_domain_seqs,
                     mibig_set_indices, mibig_set, rundata_networks_per_run,
                     html_subs_per_run, mix=False):
    """Performs pairwise comparison between BGCs. By default, this only compares BGCs from the
    input set. With --mix enabled, this also creates a mix class containing all BGCs for an
    all-vs-all comparison

    Inputs:
        run: run details for this execution of BiG-SCAPE
        bgc_collection: BgcCollection object containing all BGCs collected in this run for
            comparison
        aligned_domain_seqs: list of aligned domain sequences from hmm.read_aligned_files
        migib_set_indices: list of mibig set BGC indices used to discern mibig BGCs from
            input BGCs
        mibig_set: set of paths pointing to mibig gbk files
        rundata_networks_per_run: TODO
        html_subs_per_run: TODO:
        mix: boolean indicating whether to use a mix class. Default: False
    """
    logging.info(" Working for each BGC class")

    # we have to find the idx of query_bgc
    if run.directories.has_query_bgc:
        if run.directories.query_bgc_name in bgc_collection.bgc_name_set:
            query_bgc_idx = bgc_collection.bgc_name_tuple.index(run.directories.query_bgc_name)
        else:
            logging.error("Error finding the index of Query BGC")
            sys.exit(1)

    # we will need this for the two distance matrix generation calls

    # cluster_cache = get_cluster_cache_async(run, cluster_names, domain_list, domain_count_per_gene, corebiosynthetic_pos, bgc_gene_orientation, domain_name_info, bgc_info)
    # create working set with indices of valid clusters
    bgc_classes = create_working_set(run, bgc_collection, mix)

    if mix:
        logging.info("  %s (%d BGCs)", "Mix", len(bgc_classes["mix"]))

        # create output directory
        create_directory(os.path.join(run.directories.network, "mix"), "  Mix", False)
    else:
        # Preparing gene cluster classes
        logging.info("  Sorting the input BGCs")

        class_names_len = len(run.distance.bgc_class_names)
        bgc_class_name_2_index = dict(zip(run.distance.bgc_class_names, range(class_names_len)))

    # only make folders for the bgc_classes that are found
    for bgc_class in bgc_classes:
        if run.directories.has_query_bgc:
            # not interested in this class if our Query BGC is not here...
            if query_bgc_idx not in bgc_classes[bgc_class]:
                continue

        logging.info("  %s (%d BGCs)", bgc_class, len(bgc_classes[bgc_class]))
        if run.mibig.use_mibig:
            if len(set(bgc_classes[bgc_class]) & mibig_set_indices) == len(bgc_classes[bgc_class]):
                logging.info(" - All clusters in this class are MIBiG clusters -")
                logging.info("  If you'd like to analyze MIBiG clusters, turn off the --mibig option")
                logging.info("  and point --inputdir to the Annotated_MIBiG_reference folder")
                continue

        # create output directory
        create_directory(os.path.join(run.directories.network, bgc_class), "  All - " + bgc_class, False)

        # Create an additional file with the final list of all clusters in the class
        logging.info("   Writing annotation files")
        network_annotation_path = os.path.join(run.directories.network, bgc_class, "Network_Annotations_" + bgc_class + ".tsv")
        with open(network_annotation_path, "w") as network_annotation_file:
            network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
            for idx in bgc_classes[bgc_class]:
                bgc = bgc_collection.bgc_name_tuple[idx]
                accession_id = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.accession_id
                description = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.description
                product = bgc_collection.bgc_collection_dict[bgc].bgc_info.product
                organism = bgc_collection.bgc_collection_dict[bgc].bgc_info.organism
                taxonomy = bgc_collection.bgc_collection_dict[bgc].bgc_info.taxonomy
                network_annotation_file.write("\t".join([bgc, accession_id, description, product, sort_bgc(product), organism, taxonomy]) + "\n")

            
        if run.directories.has_query_bgc:
            pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx], bgc_classes[bgc_class])])
        else:
        # convert into a set of ordered tuples
            pairs = set([tuple(sorted(combo)) for combo in combinations(bgc_classes[bgc_class], 2)])

        if mix:
            pairs = [(x, y, -1) for (x, y) in pairs]
        else:
            pairs = [(x, y, bgc_class_name_2_index[bgc_class]) for (x, y) in pairs]

        bgc_name_id_dict = {name: idx for idx, name in enumerate(bgc_collection.bgc_name_tuple)}

        pair_count = len(pairs)
        filtered_pairs_jaccard = 0
        filtered_pairs_features = 0
        
        network_matrix = []

        # get jaccard treshold from options
        jaccard_threshold = None
        if run.options.jaccard_filter:
            jaccard_threshold = run.options.jaccard_threshold
            logging.info("    Using jaccard treshold filtering: %f", jaccard_threshold)

            remaining_pairs = 0
            remaining_pair_list = []
            for bgc_a, bgc_b, group in pairs:
                bgc_name_a = bgc_collection.bgc_name_tuple[bgc_a]
                bgc_name_b = bgc_collection.bgc_name_tuple[bgc_b]
                bgc_info_a = bgc_collection.bgc_collection_dict[bgc_name_a]
                bgc_info_b = bgc_collection.bgc_collection_dict[bgc_name_b]

                intersect = bgc_info_a.ordered_domain_set & bgc_info_b.ordered_domain_set
                overlap = bgc_info_a.ordered_domain_set | bgc_info_b.ordered_domain_set

                jaccard_idx = calc_jaccard(intersect, overlap)
                
                if jaccard_idx < jaccard_threshold:
                    network_matrix.append(generate_unrelated_row(bgc_a, bgc_b))
                    filtered_pairs_jaccard += 1
                else:
                    remaining_pair_list.append([bgc_a, bgc_b, group])
                    remaining_pairs += 1
            logging.info(
                "    %d/%d pairs with jaccard < %f filtered out",
                filtered_pairs_jaccard,
                pair_count,
                jaccard_threshold
            )
            pairs = remaining_pair_list

        # cosine distance filtering from features
        if run.options.feature_filter:
            logging.info("    Generating a list of skippable pairs using numerical features")
            logging.info("    Loading stored info from database")

            bgc_names = get_bgc_names(database)
            hmm_ids = get_hmm_ids(database)

            feature_matrix = pd.DataFrame(
                np.zeros((len(bgc_names), len(hmm_ids)), dtype=np.uint8),
                columns=hmm_ids
            )

            bgc_hmm_features = get_features(database)

            # keep track of which bgcs have features for later
            bgcs_features = set()
            # fetch feature values from db
            for bgc_name, hmm_id, value in bgc_hmm_features:
                bgc_id = bgc_name_id_dict[bgc_name]
                feature_matrix.loc[bgc_id, hmm_id] = value
                bgcs_features.add(bgc_id)

            logging.info("    Calculating cosine distances")
            cosine_dists = get_corr_cosine_dists(
                run,
                pairs,
                feature_matrix
            )
            logging.info("    Filtering pairs")
            remaining_pairs = 0
            remaining_pair_list = []
            for distance in cosine_dists:
                bgc_a_id = distance[0]
                bgc_b_id = distance[1]
                if distance[3] < run.options.feature_threshold:
                    group = distance[2]
                    remaining_pair_list.append([bgc_a_id, bgc_b_id, group])
                    remaining_pairs += 1
                else:
                    network_matrix.append(generate_unrelated_row(bgc_a_id, bgc_b_id))
                    filtered_pairs_features += 1

            # re-add anything that was lost and should be included
            for pair in pairs:
                in_features = pair[0] in bgcs_features or pair[1] in bgcs_features
                if not in_features:
                    remaining_pair_list.append(pair)
                    remaining_pairs += 1

            logging.info(
                "    %d/%d pairs with distance > %f in cosine distances filtered out",
                filtered_pairs_features,
                pair_count - filtered_pairs_jaccard,
                run.options.feature_threshold
            )
            pairs = remaining_pair_list

        if run.options.jaccard_filter and run.options.feature_filter:
            logging.info(
                "    %d/%d total pairs filtered out",
                filtered_pairs_jaccard + filtered_pairs_features,
                pair_count
            )

        # generate network matrix
        network_matrix.extend(gen_dist_matrix_async(
            run,
            database,
            pairs,
            bgc_collection,
            aligned_domain_seqs,
            jaccard_threshold
        ))

        pairs.clear()

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
                pairs = [(x, y, -1) for (x, y) in pairs]
            else:
                pairs = [(x, y, bgc_class_name_2_index[bgc_class]) for (x, y) in pairs]

            network_matrix_new_set = gen_dist_matrix_async(run, database, pairs,
                                                           bgc_collection, aligned_domain_seqs,
                                                           jaccard_threshold)
            pairs.clear()

            # Update the network matrix (QBGC-vs-all) with the distances of
            # QBGC's GCF
            network_matrix.extend(network_matrix_new_set)

            # Update actual list of BGCs that we'll use
            bgc_classes[bgc_class] = new_set
            bgc_classes[bgc_class].extend([query_bgc_idx])
            bgc_classes[bgc_class].sort()

            # Create an additional file with the list of all clusters in the class + other info
            # This version of the file only has information on the BGCs connected to Query BGC
            logging.info("   Writing annotation file (Query BGC)")
            network_annotation_path = os.path.join(run.directories.network, bgc_class, "Network_Annotations_" + bgc_class + "_QueryBGC.tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in bgc_classes[bgc_class]:
                    bgc = bgc_collection.bgc_name_tuple[idx]
                    accession_id = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.accession_id
                    description = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.description
                    product = bgc_collection.bgc_collection_dict[bgc].bgc_info.product
                    organism = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.organism
                    taxonomy = product = bgc_collection.bgc_collection_dict[bgc].bgc_info.taxonomy
                    network_annotation_file.write("\t".join([bgc, accession_id, description, product, sort_bgc(product), organism, taxonomy]) + "\n")
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

            logging.info("   Removing %d non-relevant MIBiG BGCs", len(mibig_set_del))
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

        logging.info("   Writing output files")
        path_base = os.path.join(run.directories.network, bgc_class)
        cutoffs_filenames = get_output_cutoffs_filenames(run, path_base, bgc_class)
        write_distance_matrix(network_matrix, cutoffs_filenames, run.options.include_singletons, bgc_collection)

        logging.info("  Calling Gene Cluster Families")
        reduced_network, pos_alignments = reduce_network(network_matrix)

        family_data = cluster_json_batch(database, bgc_classes[bgc_class], path_base, bgc_class,
            reduced_network, pos_alignments, bgc_collection,
            mibig_set, run.directories.pfd, run.directories.bgc_fasta,
            aligned_domain_seqs,
            cutoffs=run.cluster.cutoff_list, cluster_clans=run.options.clans, clan_cutoff=run.options.clan_cutoff,
            html_folder=run.directories.network_html)
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
