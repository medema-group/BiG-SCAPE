"""Contains methods to run the legacy comparison workflow on a bin of BGC pairs

This is a an alternate implementation that minimizes memory overhead and is generally
cleaner. TODO: profile this and compare performance with other implementation. this may
actually be faster
"""

# from python
import logging
from concurrent.futures import ProcessPoolExecutor, wait
from typing import Callable

# from other modules
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy
from src.network import BSNetwork

# from this module
from .binning import BGCBin, BGCPair
from .legacy_bins import LEGACY_BINS
from .legacy_extend import (
    legacy_needs_expand_pair,
    expand_glocal,
    check_expand,
    reset_expansion,
)
from .legacy_lcs import legacy_find_cds_lcs


def create_bin_network_edges_alt(
    bin: BGCBin, network: BSNetwork, alignment_mode: str, cores: int, callback: Callable
):  # pragma no cover
    pairs_todo = bin.num_pairs()
    logging.info("Performing distance calculation for %d pairs", pairs_todo)

    # prepare a process pool
    logging.info("Using %d cores", cores)
    with ProcessPoolExecutor(cores) as executor:
        done_pairs = 0

        running_tasks = {}

        pair_generator = bin.pairs()

        for i in range(cores * 3):
            if next_pair := next(pair_generator, None):
                startup_task = executor.submit(
                    calculate_scores_pair, (next_pair, alignment_mode, bin.label)
                )
                running_tasks[startup_task] = next_pair

        while True:
            done, not_done = wait(running_tasks, None, "FIRST_COMPLETED")

            # first quickly start new tasks
            for done_task in done:
                next_pair = next(pair_generator, None)

                if next_pair:
                    new_task = executor.submit(
                        calculate_scores_pair, (next_pair, alignment_mode, bin.label)
                    )
                    running_tasks[new_task] = next_pair

            # second loop to store results
            for done_task in done:
                task_pair = running_tasks[done_task]

                exception = done_task.exception()
                if exception:
                    raise exception

                dist, jc, ai, dss = done_task.result()

                del running_tasks[done_task]

                network.add_edge_pair(task_pair, dist=dist, jc=jc, ai=ai, dss=dss)

                done_pairs += 1
                callback(done_pairs)

            if len(running_tasks) == 0:
                break


def do_lcs_pair(pair: BGCPair, alignment_mode) -> bool:  # pragma no cover
    (
        a_start,
        a_stop,
        b_start,
        b_stop,
        reverse,
    ) = legacy_find_cds_lcs(
        pair.region_a.get_cds_with_domains(), pair.region_b.get_cds_with_domains()
    )

    # set the comparable region
    pair.comparable_region.a_start = a_start
    pair.comparable_region.a_stop = a_stop
    pair.comparable_region.b_start = b_start
    pair.comparable_region.b_stop = b_stop
    pair.comparable_region.reverse = reverse

    if legacy_needs_expand_pair(pair, alignment_mode):
        return True

    reset_expansion(pair.comparable_region)
    return False


def expand_pair(pair: BGCPair) -> bool:
    expand_glocal(pair.comparable_region)

    if not check_expand(pair.comparable_region):
        reset_expansion(pair.comparable_region)
        return True

    jc = calc_jaccard_pair(pair)
    if jc == 0.0:
        return False

    return True


def calculate_scores_pair(
    data: tuple[BGCPair, str, str]
) -> tuple[float, float, float, float]:  # pragma no cover
    pair, alignment_mode, bin_label = data

    jc = calc_jaccard_pair(pair)

    if jc == 0.0:
        return 1.0, 0.0, 0.0, 0.0

    # in the form [bool, Pair]. true bools means they need expansion, false they don't
    needs_expand = do_lcs_pair(pair, alignment_mode)

    still_related = True
    if needs_expand:
        still_related = expand_pair(pair)

    if not still_related:
        return 1.0, 0.0, 0.0, 0.0

    bin_weights = LEGACY_BINS[bin_label]["weights"]
    jc_weight, ai_weight, dss_weight, anchor_boost = bin_weights

    jaccard = calc_jaccard_pair(pair)

    adjacency = calc_ai_pair(pair)
    # mix anchor boost = 2.0
    dss = calc_dss_pair_legacy(pair, anchor_boost=anchor_boost)

    similarity = jaccard * jc_weight + adjacency * ai_weight + dss * dss_weight
    distance = 1 - similarity

    return distance, jaccard, adjacency, dss
