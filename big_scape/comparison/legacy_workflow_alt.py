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
from typing import Generator, Callable, Optional, TypeVar

# from other modules
from big_scape.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy

import big_scape.enums as bs_enums

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

T = TypeVar("T")


def batch_generator(generator: Generator[T, None, None], batch_size: int) -> list[T]:
    """Generate batches of items from a generator

    Args:
        generator: generator to generate batches from
        batch_size: int: size of the batches to generate

    Returns:
        list[Any]: batch of items
    """
    batch: list[T] = []
    for _ in range(batch_size):
        item = next(generator, None)
        if item is None:
            return batch
        batch.append(item)

    return batch


def generate_edges(
    pair_generator: RecordPairGenerator,
    alignment_mode: bs_enums.ALIGNMENT_MODE,
    cores: int,
    callback: Optional[Callable] = None,
    batch_size=100,
):  # pragma no cover
    """Generate edges for pairs generated using a pair generator

    Args:
        pair_generator (RecordPairGenerator): generator for pairs
        alignment_mode (str): alignment mode
        cores (int): number of cores to use
        callback (Optional[Callable]): callback to call when a batch is done
        batch_size (int): batch size to use

    yields:
        tuple[str, str, float, float, float, float]: tuple of (a_id, b_id, distance,
        jaccard, adjacency, dss)
    """
    # prepare a process pool
    logging.info("Using %d cores", cores)

    pairs = pair_generator.generate_pairs()

    batch_size = min(pair_generator.num_pairs(), batch_size)

    logging.info("Using batch size: %d", batch_size)

    with ProcessPoolExecutor(cores) as executor:
        done_pairs = 0

        running_tasks = {}

        for _ in range(cores):
            batch = batch_generator(pairs, batch_size)

            if len(batch) == 0:
                break

            new_task = executor.submit(
                calculate_scores_pair, (batch, alignment_mode, pair_generator.label)
            )
            running_tasks[new_task] = batch

        while True:
            done, _ = wait(running_tasks, None, "FIRST_COMPLETED")

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
                task_batch: list[RecordPair] = running_tasks[done_task]

                exception = done_task.exception()
                if exception:
                    raise exception

                results = done_task.result()

                del running_tasks[done_task]

                if len(results) != len(task_batch):
                    raise ValueError("Mismatch between task length and result length")

                for idx, pair in enumerate(task_batch):
                    distance, jaccard, adjacency, dss = results[idx]
                    yield (
                        pair.region_a._db_id,
                        pair.region_b._db_id,
                        distance,
                        jaccard,
                        adjacency,
                        dss,
                    )

                done_pairs += len(results)
                if callback:
                    callback(done_pairs)

            if len(running_tasks) == 0:
                break


def do_lcs_pair(
    pair: RecordPair, alignment_mode: bs_enums.ALIGNMENT_MODE
) -> bool:  # pragma no cover
    """Find the longest common subsequence of protein domains between two regions

    Args:
        pair (RecordPair): pair to find the lcs for

    Returns:
        bool: True if the pair needs expansion, False if it does not
    """
    a_start, a_stop, b_start, b_stop, reverse = legacy_find_cds_lcs(
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
    """Expand the pair and calculate the jaccard index

    Args:
        pair (RecordPair): pair to expand

    Returns:
        float: jaccard index
    """
    expand_glocal(pair.comparable_region)

    if not check_expand(pair.comparable_region):
        reset_expansion(pair.comparable_region)
        jc = calc_jaccard_pair(pair)
        return jc

    jc = calc_jaccard_pair(pair)

    return jc


def calculate_scores_pair(
    data: tuple[list[RecordPair], bs_enums.ALIGNMENT_MODE, str]
) -> list[tuple[float, float, float, float]]:  # pragma no cover
    """Calculate the scores for a list of pairs

    Args:
        data (tuple[list[RecordPair], str, str]): list of pairs, alignment mode, bin
        label

    Returns:
        list[tuple[float, float, float, float]]: list of scores for each pair in the
        order as the input data list
    """
    pairs, alignment_mode, bin_label = data

    results = []

    for pair in pairs:
        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            results.append((1.0, 0.0, 0.0, 0.0))
            continue

        # in the form [bool, Pair]. true bools means they need expansion, false they don't
        needs_expand = do_lcs_pair(pair, alignment_mode)

        if needs_expand:
            jaccard = expand_pair(pair)

        if jaccard == 0.0:
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