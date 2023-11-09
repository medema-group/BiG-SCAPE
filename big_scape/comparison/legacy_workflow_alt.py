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
from math import ceil


# from other modules
from big_scape.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy
import big_scape.enums as bs_enums

# from this module
from .binning import RecordPairGenerator, RecordPair
from .binning import LEGACY_WEIGHTS
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


def get_batch_size(cores: int, desired_batch_size: int, num_items: int):
    """Get the batch size to use for a given number of cores

    The idea is to maximize the number of pairs for each core while keeping the batch
    size reasonable

    Args:
        cores (int): number of cores to use
        desired_batch_size (int): desired batch size
        num_items (int): number of items to generate batches for

    Returns:
        int: batch size to use
    """

    if num_items < cores:
        return num_items

    if num_items < desired_batch_size:
        return num_items

    if cores * desired_batch_size > num_items:
        return ceil(num_items / cores)

    return desired_batch_size


def generate_edges(
    pair_generator: RecordPairGenerator,
    alignment_mode: bs_enums.ALIGNMENT_MODE,
    cores: int,
    callback: Optional[Callable] = None,
    batch_size=None,
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
    logging.debug("Using %d cores", cores)

    pairs = pair_generator.generate_pairs()

    num_pairs = pair_generator.num_pairs()

    if batch_size is None:
        batch_size = get_batch_size(cores, 10000, num_pairs)
        logging.debug("Using automatic batch size: %d", batch_size)
    else:
        batch_size = min(num_pairs, batch_size)
        logging.debug("Using batch size: %d", batch_size)

    with ProcessPoolExecutor(cores) as executor:
        done_pairs = 0

        running_tasks = {}

        for _ in range(cores):
            batch = batch_generator(pairs, batch_size)

            if len(batch) == 0:
                break

            new_task = executor.submit(
                calculate_scores_pair,
                (
                    batch,
                    alignment_mode,
                    pair_generator.edge_param_id,
                    pair_generator.weights,
                ),
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
                    calculate_scores_pair,
                    (
                        batch,
                        alignment_mode,
                        pair_generator.edge_param_id,
                        pair_generator.weights,
                    ),
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

                yield from results

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
        pair.record_a.get_cds_with_domains(), pair.record_b.get_cds_with_domains()
    )

    # set the comparable region
    pair.comparable_region.lcs_a_start = a_start
    pair.comparable_region.lcs_a_stop = a_stop
    pair.comparable_region.lcs_b_start = b_start
    pair.comparable_region.lcs_b_stop = b_stop
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

    pair.comparable_region.alignment_mode = bs_enums.ALIGNMENT_MODE.GLOCAL
    jc = calc_jaccard_pair(pair)

    return jc


def calculate_scores_pair(
    data: tuple[list[RecordPair], bs_enums.ALIGNMENT_MODE, int, str]
) -> list[
    tuple[
        Optional[int],
        Optional[int],
        float,
        float,
        float,
        float,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        bool,
    ]
]:  # pragma no cover
    """Calculate the scores for a list of pairs

    Args:
        data (tuple[list[RecordPair], str, str]): list of pairs, alignment mode, bin
        label

    Returns:
        list[tuple[int, int, float, float, float, float, int, int, int, int, int, int,
        int, int, bool, str,]]: list of scores for each pair in the
        order as the input data list, including lcs and extension coordinates
    """
    pairs, alignment_mode, edge_param_id, weights_label = data

    results = []

    # TODO: this fails since DB getting accessed from child processes
    # seems to be a problem with the DB connection (for mac?)
    # weights_label = bs_comparison.get_edge_weight(edge_param_id)

    for pair in pairs:
        if pair.record_a.parent_gbk == pair.record_b.parent_gbk:
            results.append(
                (
                    pair.record_a._db_id,
                    pair.record_b._db_id,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    edge_param_id,
                    pair.comparable_region.lcs_a_start,
                    pair.comparable_region.lcs_a_stop,
                    pair.comparable_region.lcs_b_start,
                    pair.comparable_region.lcs_b_stop,
                    pair.comparable_region.a_start,
                    pair.comparable_region.a_stop,
                    pair.comparable_region.b_start,
                    pair.comparable_region.b_stop,
                    pair.comparable_region.reverse,
                )
            )
            continue

        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            results.append(
                (
                    pair.record_a._db_id,
                    pair.record_b._db_id,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    edge_param_id,
                    pair.comparable_region.lcs_a_start,
                    pair.comparable_region.lcs_a_stop,
                    pair.comparable_region.lcs_b_start,
                    pair.comparable_region.lcs_b_stop,
                    pair.comparable_region.a_start,
                    pair.comparable_region.a_stop,
                    pair.comparable_region.b_start,
                    pair.comparable_region.b_stop,
                    pair.comparable_region.reverse,
                )
            )
            continue

        # in the form [bool, Pair]. true bools means they need expansion, false they don't
        needs_expand = do_lcs_pair(pair, alignment_mode)

        if needs_expand:
            jaccard = expand_pair(pair)

        if jaccard == 0.0:
            results.append(
                (
                    pair.record_a._db_id,
                    pair.record_b._db_id,
                    1.0,
                    0.0,
                    0.0,
                    0.0,
                    edge_param_id,
                    pair.comparable_region.lcs_a_start,
                    pair.comparable_region.lcs_a_stop,
                    pair.comparable_region.lcs_b_start,
                    pair.comparable_region.lcs_b_stop,
                    pair.comparable_region.a_start,
                    pair.comparable_region.a_stop,
                    pair.comparable_region.b_start,
                    pair.comparable_region.b_stop,
                    pair.comparable_region.reverse,
                )
            )
            continue

        if weights_label not in LEGACY_WEIGHTS:
            bin_weights = LEGACY_WEIGHTS["mix"]["weights"]
        else:
            bin_weights = LEGACY_WEIGHTS[weights_label]["weights"]
        jc_weight, ai_weight, dss_weight, anchor_boost = bin_weights

        jaccard = calc_jaccard_pair(pair)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=anchor_boost)

        similarity = jaccard * jc_weight + adjacency * ai_weight + dss * dss_weight
        distance = 1 - similarity

        results.append(
            (
                pair.record_a._db_id,
                pair.record_b._db_id,
                distance,
                jaccard,
                adjacency,
                dss,
                edge_param_id,
                pair.comparable_region.lcs_a_start,
                pair.comparable_region.lcs_a_stop,
                pair.comparable_region.lcs_b_start,
                pair.comparable_region.lcs_b_stop,
                pair.comparable_region.a_start,
                pair.comparable_region.a_stop,
                pair.comparable_region.b_start,
                pair.comparable_region.b_stop,
                pair.comparable_region.reverse,
            )
        )

    return results
