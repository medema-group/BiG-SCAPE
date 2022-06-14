"""Module containing code for orchestrating the multithreader distance
calculation

Authors: Jorge Navarro, Arjan Draisma
"""
import logging
import math
from array import array

from multiprocessing import Pipe, Process
from multiprocessing.connection import Connection, wait

from src.big_scape.bgc_dom_info import BgcDomainInfo
from src.big_scape.bgc_collection import BgcCollection
from src.big_scape.bgc_info import BgcInfo
from src.big_scape.scores import calc_adj_idx, calc_distance, calc_dss, calc_jaccard, gen_unrelated_pair_distance, process_orientation

def gen_dist_matrix_worker(
    connection: Connection,
    run,
    database,
    bgc_collection,
    aligned_domain_sequences,
    jaccard_threshold
):
    """Worker method for threads that process distance calculation. Takes bgc
    pairs from an input queue and generates distances between the bgc pairs 
    which it puts back into the output queue.

    Inputs in the queue are expected to be a tuple of the format (pair,
    skip_set), where:
    - pair is a tuple containing two bgc indices

    This thread continues to run until an explicit exit signal is received.

    Inputs:
        connection: worker end of the pipe used for I/O
        run: run details for this execution of BiG-SCAPE
        database: database object for this run of BiG-SCAPE
        bgc_collection: collection of BGCs
        aligned_domain_sequences: list of aligned domain sequences from
                                  hmm.read_aligned_files
        jaccard_threshold: TODO: remove
    """
    while True:
        input_task_batch = connection.recv()
        results = []
        for input_task in input_task_batch:
            if input_task is None:
                if len(results) > 0:
                    connection.send(results)
                connection.close()
                return
            # logging.info("launching task on pair %s, %s", pair[0], pair[1])

            result = generate_dist_matrix(
                input_task,
                database,
                run,
                bgc_collection,
                aligned_domain_sequences,
                jaccard_threshold = jaccard_threshold
            )
            results.append(result)
        connection.send(results)


# @timeit
def gen_dist_matrix_async(
    run,
    database,
    cluster_pairs,
    bgc_collection: BgcCollection,
    aligned_domain_sequences,
    jaccard_threshold
):
    """Distributes the distance calculation part
    cluster_pairs is a list of triads (cluster1_index, cluster2_index, 
    BGC class)
    """

    num_processes = run.options.cores

    connections = list()

    num_tasks = len(cluster_pairs)

    cluster_idx = 0

    # how many tasks to send at once to a worker
    max_batch = 200

    processes = []
    for thread_num in range(num_processes):
        # if there are fewer tasks than there are threads, don't start a new
        # thread. this is only a concern for < cores tasks
        if cluster_idx >= num_tasks:
            continue
        # generate thread name
        thread_name = f"distance_thread_{thread_num}"
        logging.debug("Starting %s", thread_name)

        # generate connection
        main_connection, worker_connection = Pipe(True)
        connections.append(main_connection)

        # create and start new process
        process_args = (
            worker_connection,
            run,
            database,
            bgc_collection,
            aligned_domain_sequences,
            jaccard_threshold
        )
        new_process = Process(
            target=gen_dist_matrix_worker,
            args=process_args,
            name=thread_name
        )
        processes.append(new_process)
        new_process.start()

        # send starting data
        tasks_left = num_tasks - cluster_idx
        if tasks_left < max_batch:
            num_batch = tasks_left
            batch = cluster_pairs[cluster_idx:cluster_idx + num_batch]
            batch.append(None)
        else:
            num_batch = max_batch
            batch = cluster_pairs[cluster_idx:cluster_idx + num_batch]

        main_connection.send(batch)
        cluster_idx += num_batch


    network_matrix = []

    while True:
        all_tasks_put = cluster_idx == num_tasks
        all_tasks_done = len(network_matrix) == num_tasks

        if all_tasks_put and all_tasks_done:
            break

        available_connections = wait(connections, timeout=1)
        for connection in available_connections:
            connection: Connection

            # receive results
            network_rows = connection.recv()
            
            # close connection if no more tasks need to be added
            if all_tasks_put:
                connection.send([None])
                connection.close()
                connections.remove(connection)
            else:
                # add a new batch of tasks
                tasks_left = num_tasks - cluster_idx
                if tasks_left < max_batch:
                    num_batch = tasks_left
                else:
                    num_batch = max_batch

                batch = cluster_pairs[cluster_idx:cluster_idx + num_batch]

                # send batch
                connection.send(batch)

                # increment index
                cluster_idx += num_batch
                all_tasks_put = cluster_idx == num_tasks

            # add row to matrix
            for network_row in network_rows:
                network_matrix.append(network_row)

                num_tasks_done = len(network_matrix)
                all_tasks_done = len(network_matrix) == num_tasks

                # print progress every 10%
                if num_tasks_done % math.ceil(num_tasks / 10) == 0:
                    percent_done = num_tasks_done / num_tasks * 100
                    logging.info(
                        "    %d%% (%d/%d)",
                        percent_done,
                        num_tasks_done,
                        num_tasks
                    )


    for process in processes:
        process.join()
        thread_name = process.name
        logging.debug("Thread %s stopped", thread_name)

    return network_matrix


def calc_ai_pair(
    cluster_a: BgcInfo,
    cluster_b: BgcInfo,
    pair_dom_info: BgcDomainInfo
):
    """Calculate the adjacency index of a pair of BGCs

    Inputs:
        cluster_a, cluster_b: BgcInfo objects of the clusters to be compared
        pair_dom_info: BgcDomainInfo of this pair
    """
    a_dom_list = cluster_a.ordered_domain_list
    b_dom_list = cluster_b.ordered_domain_list

    a_dom_start = pair_dom_info.a_dom_start
    a_dom_end = pair_dom_info.a_dom_end

    b_dom_start = pair_dom_info.b_dom_start
    b_dom_end = pair_dom_info.b_dom_end

    return calc_adj_idx(
        a_dom_list,
        b_dom_list,
        a_dom_start,
        a_dom_end,
        b_dom_start,
        b_dom_end
    )

def generate_dist_matrix(
    parms,
    database,
    run,
    bgc_collection: BgcCollection,
    aligned_domain_sequences,
    jaccard_threshold = 0.0
):
    """Unpack data to actually launch cluster_distance for one pair of BGCs

    Inputs:
        parms: a tuple of the form (bgc_a_idx, bgc_b_idx, class_idx) for this
            comparison
        bgc_collection: BgcCollection object containing bgc info objects and
            necessary data for
            distance calculation
        aligned_domain_sequences: list of aligned domain sequences from
            hmm.read_aligned_files
    """

    cluster_1_idx, cluster_2_idx, bgc_class_idx = parms

    bgc_class = run.distance.bgc_class_names[bgc_class_idx]

    cluster_name_a = bgc_collection.bgc_name_tuple[cluster_1_idx]
    cluster_name_b = bgc_collection.bgc_name_tuple[cluster_2_idx]

    cluster_a = bgc_collection.bgc_collection_dict[cluster_name_a]
    cluster_b = bgc_collection.bgc_collection_dict[cluster_name_b]

    # this really shouldn't happen if we've filtered domain-less gene clusters
    # already
    if len(cluster_a.ordered_domain_list) == 0 or len(cluster_b.ordered_domain_list) == 0:
        logging.warning("   Regarding distance between clusters %s and %s:", cluster_name_a, cluster_name_b)
        if len(cluster_a.ordered_domain_list) == 0 and len(cluster_b.ordered_domain_list) == 0:
            logging.warning("   None have identified domains. Distance cannot be calculated")
        elif (cluster_a.ordered_domain_list) == 0:
            logging.warning("   Cluster %s has no identified domains. Distance set to 1", cluster_name_a)
        else:
            logging.warning("   Cluster %s has no identified domains. Distance set to 1", cluster_name_b)

        # last two values (S, Sa) should really be zero but this could give rise to errors when parsing
        # the network file (unless we catched the case S = Sa = 0

        # cluster1Idx, cluster2Idx, distance, jaccard, DSS, AI, rDSSNa, rDSSa,
        #   S, Sa, lcsStartA, lcsStartB
        return array('f', [cluster_1_idx, cluster_2_idx, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0])

    # unpack weights
    weights = run.distance.bgc_class_weight[bgc_class]
    anchor_boost = weights[3]

    # initialize domain specific info
    # this contains the domain slice information, which may change if expansion is needed
    pair_dom_info = BgcDomainInfo(cluster_a, cluster_b)

    # Detect totally unrelated pairs from the beginning
    # lack of intersect
    no_intersect = len(pair_dom_info.intersect) == 0

    if no_intersect:
        score_data = gen_unrelated_pair_distance(run, cluster_a, cluster_b)
        jaccard_index, dss, adjacency_index, dss_non_anchor, dss_anchor, num_non_anchor_domains, num_anchor_domains, slice_start_a, slice_start_b, slice_length_a, rev = score_data
    else:
        slice_data = process_orientation(cluster_a, cluster_b)

        slice_start_a, slice_start_b, slice_length_a, slice_length_b, use_b_string, slice_reverse = slice_data

        rev = 0.0
        if slice_reverse:
            rev = 1.0

        pair_dom_info.expand_score(run, cluster_a, cluster_b, slice_data)


        # JACCARD INDEX
        union = pair_dom_info.a_dom_set | pair_dom_info.b_dom_set
        jaccard_index = calc_jaccard(pair_dom_info.intersect, union)

        dss_data = calc_dss(
            run,
            database,
            cluster_a,
            cluster_b,
            aligned_domain_sequences,
            anchor_boost,
            pair_dom_info
        )
        # unpack variables
        dss, dss_non_anchor, dss_anchor, num_non_anchor_domains, num_anchor_domains = dss_data

        adjacency_index = calc_ai_pair(cluster_a, cluster_b, pair_dom_info)


    dist = calc_distance(weights, jaccard_index, dss, adjacency_index, cluster_a.name, cluster_b.name)


    network_row = array('f', [cluster_1_idx, cluster_2_idx, dist, (1-dist)**2, jaccard_index,
                              dss, adjacency_index, dss_non_anchor, dss_anchor, num_non_anchor_domains, num_anchor_domains, slice_start_a, slice_start_b,
                              slice_length_a, rev])
    return network_row



def write_distance_matrix(distance_matrix, cutoffs_filenames, include_singletons, bgc_collection: BgcCollection):
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
        gene_cluster_a = bgc_collection.bgc_name_tuple[int(matrix_entry[0])]
        gene_cluster_b = bgc_collection.bgc_name_tuple[int(matrix_entry[1])]
        row = [gene_cluster_a, gene_cluster_b]

        # get AntiSMASH annotations
        cluster_group_a = bgc_collection.bgc_collection_dict[gene_cluster_a].bgc_info.product
        cluster_group_b = bgc_collection.bgc_collection_dict[gene_cluster_b].bgc_info.product

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
