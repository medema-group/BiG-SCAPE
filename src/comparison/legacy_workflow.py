"""Contains methods to run the legacy comparison workflow on a bin of BGC pairs"""

# from python
import logging
from math import ceil
from multiprocessing import cpu_count
from multiprocessing.connection import Connection, wait
from typing import cast

# from other modules
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy
from src.network import BSNetwork
from src.utility import start_processes

# from this module
from .legacy_extend import (
    legacy_needs_extend,
    expand_glocal,
    check_expand,
    reset_expansion,
)
from .binning import BGCBin, BGCPair


def create_bin_network_edges(bin: BGCBin, network: BSNetwork, alignment_mode: str):
    # first step is to calculate the Jaccard of all pairs. This is pretty fast, but
    # could be optimized by multiprocessing for very large bins
    logging.info("Calculating Jaccard for %d pairs", bin.num_pairs())

    for pair in bin.pairs(legacy_sorting=True):
        # calculate jaccard for the full sets. if this is 0, there are no shared domains
        # important not to cache here otherwise we are using the full range again later
        jaccard = calc_jaccard_pair(pair, cache=False)

        if jaccard == 0.0:
            network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

    # any pair that had a jaccard of 0 are put into the network and should not be
    # processed again

    # next step is to perform LCS. We need to multiprocess this and that is a bit of a
    # hassle

    logging.info(
        "Performing LCS for %d pairs", bin.num_pairs() - network.graph.number_of_edges()
    )

    pairs_need_expand = []
    pairs_no_expand = []
    for pair in bin.pairs(legacy_sorting=True):
        if pair in network:
            continue

        logging.debug("JC: %f", jaccard)

        pair.comparable_region.find_lcs()
        pair.comparable_region.log_comparable_region("LCS")

        if legacy_needs_extend(pair, alignment_mode):
            pairs_need_expand.append(pair)
            continue

        # or it is not expanded at all and the entire region is used
        reset_expansion(pair.comparable_region)
        pairs_no_expand.append(pair)

    # those regions which need expansion are now expanded. Expansion is expensive and
    # is also done through multiprocessing
    logging.info("Expanding regions for %d pairs", len(pairs_need_expand))

    expanded_pairs = []
    for pair in pairs_need_expand:
        expand_glocal(pair.comparable_region)

        # if after expansion the region is still too small or does not contain any
        # biosynthetic genes, we reset back to the full region and add this pair
        # to the list of pairs that were not expanded
        if not check_expand(pair.comparable_region):
            reset_expansion(pair.comparable_region)
            pairs_no_expand.append(pair)
            continue

        pair.comparable_region.log_comparable_region("GLOCAL")

        jaccard = calc_jaccard_pair(pair)

        # any pair with a jaccard of 0 after expansion is also kicked out
        if jaccard == 0.0:
            network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

        expanded_pairs.append(pair)

    # from here on the only things left to be done

    logging.info(
        "Calculating score for %d pairs that were reset or did not need expansion",
        len(pairs_no_expand),
    )
    calculate_scores_multiprocess(pairs_no_expand, network)

    logging.info(
        "Calculating score for %d pairs that were expanded",
        len(expanded_pairs),
    )
    calculate_scores_multiprocess(expanded_pairs, network)


def calculate_scores_worker_method(
    pair_idx: int, pairs: list[BGCPair]
) -> tuple[int, float, float, float, float]:
    """Calculate and return the scores and distance for a pair

    Args:
        pair_idx: The index in the original list of this pair.
        pair (BGCPair): Pair of BGCs to calculate scores for

    Returns:
        tuple[float, float, float, float]: distance, jaccard, AI, DSS
    """
    pair = pairs[pair_idx]

    jaccard = calc_jaccard_pair(pair)

    adjacency = calc_ai_pair(pair)
    # mix anchor boost = 2.0
    dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

    # mix
    distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

    return pair_idx, distance, jaccard, adjacency, dss


def calculate_scores_multiprocess(
    pairs: list[BGCPair],
    network: BSNetwork,
    num_processes: int = cpu_count(),
    batch_size=None,
):
    """Calculate the scores for a list of pairs by using subprocesses

    Args:
        pairs (list[BGCPair]): list of pairs to perform score calution on
        network (BSNetwork): BSnetwork objects to add new edges to
        cpu_count (int): number of cores to use. Defaults to number of cpus available
        batch_size (int): size of the batches to send to the worker. IF set to none,
        evenly divides the task set into number of batches equal to num_processes
    """

    # prepare processes
    processes, connections = start_processes(
        num_processes, calculate_scores_worker_method, pairs
    )

    if batch_size is None:
        batch_size = ceil(len(pairs) / num_processes)

    tasks_done = 0
    pair_idx = 0

    while len(connections) > 0:
        available_connections = wait(connections)

        for connection in available_connections:
            connection = cast(Connection, connection)

            output_data = connection.recv()

            if pair_idx < len(pairs):
                sent_task_num = min(batch_size, len(pairs) - pair_idx)
                input_data = [sent_task_num]
                input_data.extend(range(pair_idx, pair_idx + sent_task_num))
                pair_idx += sent_task_num
            else:
                input_data = None

            connection.send(input_data)

            if input_data is None:
                connection.close()
                connections.remove(connection)

            if output_data is not None:
                recv_task_num = output_data[0]
                for task_output in output_data[1:]:
                    done_pair_idx, dist, jc, ai, dss = task_output

                    network.add_edge_pair(
                        pairs[done_pair_idx], dist=dist, jc=jc, ai=ai, dss=dss
                    )

                tasks_done += recv_task_num

    # just to make sure, kill any remaining processes
    for process in processes:
        process.kill()
