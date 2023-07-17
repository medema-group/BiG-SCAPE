"""Contains methods to run the legacy comparison workflow on a bin of BGC pairs"""

# from python
import logging
from multiprocessing import Pipe, Process, cpu_count
from multiprocessing.connection import Connection, wait
from typing import cast

# from other modules
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy
from src.network import BSNetwork

# from this module
from .legacy_extend import (
    legacy_needs_extend,
    expand_glocal,
    check_expand,
    reset_expansion,
)
from .binning import BGCBin, BGCPair


def create_bin_network_edges_old(bin: BGCBin, network: BSNetwork, alignment_mode: str):
    for pair in bin.pairs(legacy_sorting=True):
        # calculate jaccard for the full sets. if this is 0, there are no shared domains
        # important not to cache here otherwise we are using the full range again later
        jaccard = calc_jaccard_pair(pair, cache=False)

        if jaccard == 0.0:
            network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

        logging.debug("JC: %f", jaccard)

        # we record the LCS starts and stops here in case we need to reset them later
        pair.comparable_region.find_lcs()
        pair.comparable_region.log_comparable_region("LCS")

        if legacy_needs_extend(pair, alignment_mode):
            expand_glocal(pair.comparable_region)

            if check_expand(pair.comparable_region):
                pair.comparable_region.log_comparable_region("GLOCAL")

                jaccard = calc_jaccard_pair(pair)

                if jaccard == 0.0:
                    network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
                    continue
            else:
                reset_expansion(pair.comparable_region)
        else:
            reset_expansion(pair.comparable_region)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

        # mix
        distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

        logging.debug(
            "JC: %f, AI: %f, DSS: %f, SCORE: %f", jaccard, adjacency, dss, distance
        )

        network.add_edge_pair(pair, jc=jaccard, ai=adjacency, dss=dss, dist=distance)


def create_bin_network_edges(bin: BGCBin, network: BSNetwork, alignment_mode: str):
    logging.info("Calculating Jaccard for %d pairs", bin.num_pairs())
    for pair in bin.pairs(legacy_sorting=True):
        # calculate jaccard for the full sets. if this is 0, there are no shared domains
        # important not to cache here otherwise we are using the full range again later
        jaccard = calc_jaccard_pair(pair, cache=False)

        if jaccard == 0.0:
            network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

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

        reset_expansion(pair.comparable_region)
        pairs_no_expand.append(pair)

    logging.info("Expanding regions for %d pairs", len(pairs_need_expand))
    expanded_pairs = []
    for pair in pairs_need_expand:
        expand_glocal(pair.comparable_region)

        if not check_expand(pair.comparable_region):
            reset_expansion(pair.comparable_region)
            pairs_no_expand.append(pair)
            continue

        pair.comparable_region.log_comparable_region("GLOCAL")

        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            network.add_edge_pair(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

        expanded_pairs.append(pair)

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


def calculate_scores_multiprocess(pairs: list[BGCPair], network: BSNetwork):
    """Calculate the scores for a list of pairs by using subprocesses

    Args:
        pairs (list[BGCPair]): list of pairs to perform score calution on
        network (BSNetwork): BSnetwork objects to add new edges to
    """

    # prepare processes
    processes = []
    connections: list[Connection] = []
    for process_id in range(cpu_count()):
        main_connection, worker_connection = Pipe()
        connections.append(main_connection)

        process = Process(
            target=calculate_pair_worker, args=(process_id, worker_connection, pairs)
        )
        processes.append(process)
        process.start()

    tasks_done = 0
    pair_idx = 0

    while len(connections) > 0:
        available_connections = wait(connections)

        for connection in available_connections:
            connection = cast(Connection, connection)

            output_data = connection.recv()

            if pair_idx < len(pairs):
                input_data = pair_idx
                pair_idx += 1
            else:
                input_data = None

            connection.send(input_data)

            if input_data is None:
                connection.close()
                connections.remove(connection)

            if output_data is not None:
                done_pair_idx, dist, jc, ai, dss = output_data

                tasks_done += 1

                network.add_edge_pair(
                    pairs[done_pair_idx], dist=dist, jc=jc, ai=ai, dss=dss
                )

    # just to make sure, kill any remaining processes
    for process in processes:
        process.kill()


def calculate_pair_worker(
    process_id: int, worker_connection: Connection, pairs: list[BGCPair]
) -> None:
    """Waits for tasks to arrive from the main process and returns pair distances through a
    connection

    Args:
        process_id (int): ID of this process
        worker_connection (Connection): Worker connection to connect to the main process
        pairs (list[BGCPair]): List of pairs to refer to
    """

    worker_connection.send(None)

    while True:
        task = worker_connection.recv()

        if task is None:
            return

        pair = pairs[task]

        scores = calculate_pair_scores(pair)

        output = (task,) + scores

        worker_connection.send(output)


def calculate_pair_scores(pair: BGCPair) -> tuple[float, float, float, float]:
    """Calculate and return the scores and distance for a pair

    Args:
        pair (BGCPair): Pair of BGCs to calculate scores for

    Returns:
        tuple[float, float, float, float]: distance, jaccard, AI, DSS
    """
    jaccard = calc_jaccard_pair(pair)

    adjacency = calc_ai_pair(pair)
    # mix anchor boost = 2.0
    dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

    # mix
    distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

    return distance, jaccard, adjacency, dss
