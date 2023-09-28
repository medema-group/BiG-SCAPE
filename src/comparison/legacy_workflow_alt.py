"""Contains methods to run the legacy comparison workflow on a bin of BGC pairs

This is a an alternate implementation that minimizes memory overhead and is generally
cleaner.
TODO: profile this and compare performance with other implementation. this may
actually be faster
TODO: this is relatively simple code and worth testing, except
create_bin_network_edges_alt
TODO: docstrings, none typings
"""

# from python
import logging
from concurrent.futures import ProcessPoolExecutor, wait
from typing import Callable, Optional

# from other modules
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy

# from this module
from .binning import RecordPairGenerator, RecordPair
from .legacy_bins import LEGACY_BINS
from .legacy_extend import (
    legacy_needs_expand_pair,
    expand_glocal,
    check_expand,
    reset_expansion,
)
from .legacy_lcs import legacy_find_cds_lcs


def batch_generator(iterator, batch_size):
    batch = []
    for i in range(batch_size):
        item = next(iterator, None)
        if item is None:
            return batch
        batch.append(item)

    return batch


def generate_edges(
    pair_generator: RecordPairGenerator,
    alignment_mode: str,
    cores: int,
    callback: Optional[Callable] = None,
    batch_size=100,
):  # pragma no cover
    # prepare a process pool
    logging.info("Using %d cores", cores)

    pairs = pair_generator.generate_pairs()

    batch_size = min(pair_generator.num_pairs(), batch_size)

    logging.info("Using batch size: %d", batch_size)

    with ProcessPoolExecutor(cores) as executor:
        done_pairs = 0

        running_tasks = {}

        for i in range(cores):
            batch = batch_generator(pairs, batch_size)

            if len(batch) == 0:
                break

            new_task = executor.submit(
                calculate_scores_pair, (batch, alignment_mode, pair_generator.label)
            )
            running_tasks[new_task] = batch

        while True:
            done, not_done = wait(running_tasks, None, "FIRST_COMPLETED")

            # first quickly start new tasks
            for done_task in done:
                batch = batch_generator(pairs, batch_size)

                if len(batch) == 0:
                    continue

                new_task = executor.submit(
                    calculate_scores_pair, (batch, alignment_mode, pair_generator.label)
                )
                running_tasks[new_task] = batch

            # second loop to store results
            for done_task in done:
                task_batch = running_tasks[done_task]

                exception = done_task.exception()
                if exception:
                    raise exception

                results = done_task.result()

                del running_tasks[done_task]

                if len(results) != len(task_batch):
                    raise ValueError("Mismatch between task length and result length")

                for idx, pair in enumerate(task_batch):
                    dist, jc, ai, dss = results[idx]
                    yield (
                        pair.region_a._db_id,
                        pair.region_b._db_id,
                        dist,
                        jc,
                        ai,
                        dss,
                    )

                done_pairs += len(results)
                if callback:
                    callback(done_pairs)

            if len(running_tasks) == 0:
                break


def do_lcs_pair(pair: RecordPair, alignment_mode) -> bool:  # pragma no cover
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


def expand_pair(pair: RecordPair) -> float:
    expand_glocal(pair.comparable_region)

    if not check_expand(pair.comparable_region):
        reset_expansion(pair.comparable_region)
        jc = calc_jaccard_pair(pair)
        return jc

    jc = calc_jaccard_pair(pair)

    return jc


def calculate_scores_pair(
    data: tuple[list[RecordPair], str, str]
) -> list[tuple[float, float, float, float]]:  # pragma no cover
    pairs, alignment_mode, bin_label = data

    results = []

    for pair in pairs:
        jc = calc_jaccard_pair(pair)

        if jc == 0.0:
            results.append((1.0, 0.0, 0.0, 0.0))
            continue

        # in the form [bool, Pair]. true bools means they need expansion, false they don't
        needs_expand = do_lcs_pair(pair, alignment_mode)

        if needs_expand:
            jc = expand_pair(pair)

        if jc == 0.0:
            results.append((1.0, 0.0, 0.0, 0.0))
            continue

        if bin_label not in LEGACY_BINS:
            bin_weights = LEGACY_BINS["mix"]["weights"]
        else:
            bin_weights = LEGACY_BINS[bin_label]["weights"]
        jc_weight, ai_weight, dss_weight, anchor_boost = bin_weights

        jaccard = calc_jaccard_pair(pair)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=anchor_boost)

        similarity = jaccard * jc_weight + adjacency * ai_weight + dss * dss_weight
        distance = 1 - similarity

        results.append((distance, jaccard, adjacency, dss))

    return results
