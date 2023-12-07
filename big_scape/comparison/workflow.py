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
from concurrent.futures import ProcessPoolExecutor, Future
from threading import Event, Condition
from typing import Generator, Callable, Optional, TypeVar
from math import ceil
from .record_pair import RecordPair


# from other modules
from big_scape.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair
import big_scape.enums as bs_enums
import big_scape.genbank as bs_gbk
from big_scape.cli.config import BigscapeConfig
import big_scape.comparison as bs_comparison

# from this module
from .binning import RecordPairGenerator
from .binning import LEGACY_WEIGHTS
from .extend import extend, reset, len_check, biosynthetic_check
from .lcs import find_domain_lcs_region, find_domain_lcs_protocluster

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
    max_queue_length: int,
    callback: Optional[Callable] = None,
    batch_size=None,
):  # pragma no cover
    """Generate edges for pairs generated using a pair generator

    Args:
        pair_generator (RecordPairGenerator): generator for pairs
        alignment_mode (str): alignment mode
        cores (int): number of cores to use
        callback (Optional[Callable]): callback to call when a batch is done. this is
        called with the results of the batch as the only argument
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

    running_futures = set()
    all_done = Event()
    new_results = Condition()

    def on_complete(future: Future):
        exception = future.exception()
        if exception is not None:
            raise exception
        with new_results:
            new_results.notify_all()
            if callback:
                callback(future.result())
            running_futures.discard(future)

    with ProcessPoolExecutor(cores) as executor:
        while True:
            while len(running_futures) < max_queue_length:
                batch = batch_generator(pairs, batch_size)

                if len(batch) == 0:
                    all_done.set()
                    break

                future = executor.submit(
                    calculate_scores_pair,
                    (
                        batch,
                        alignment_mode,
                        pair_generator.edge_param_id,
                        pair_generator.weights,
                    ),
                )

                future.add_done_callback(on_complete)
                running_futures.add(future)

            with new_results:
                new_results.wait()

            if all_done.is_set():
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

    if isinstance(pair.record_a, bs_gbk.ProtoCluster) and isinstance(
        pair.record_b, bs_gbk.ProtoCluster
    ):
        logging.debug("Using protocluster lcs")
        (
            a_start,
            a_stop,
            b_start,
            b_stop,
            a_cds_start,
            a_cds_stop,
            b_cds_start,
            b_cds_stop,
            reverse,
        ) = find_domain_lcs_protocluster(pair)

    else:
        logging.debug("Using region lcs")
        (
            a_start,
            a_stop,
            b_start,
            b_stop,
            a_cds_start,
            a_cds_stop,
            b_cds_start,
            b_cds_stop,
            reverse,
        ) = find_domain_lcs_region(pair)

    logging.debug("before lcs:")
    logging.debug(pair.comparable_region)

    # set the lcs and comparable region based on domains and on cds
    pair.comparable_region.lcs_domain_a_start = a_start
    pair.comparable_region.lcs_domain_a_stop = a_stop
    pair.comparable_region.lcs_domain_b_start = b_start
    pair.comparable_region.lcs_domain_b_stop = b_stop
    pair.comparable_region.lcs_a_start = a_cds_start
    pair.comparable_region.lcs_a_stop = a_cds_stop
    pair.comparable_region.lcs_b_start = b_cds_start
    pair.comparable_region.lcs_b_stop = b_cds_stop
    pair.comparable_region.domain_a_start = a_start
    pair.comparable_region.domain_a_stop = a_stop
    pair.comparable_region.domain_b_start = b_start
    pair.comparable_region.domain_b_stop = b_stop
    pair.comparable_region.a_start = a_cds_start
    pair.comparable_region.a_stop = a_cds_stop
    pair.comparable_region.b_start = b_cds_start
    pair.comparable_region.b_stop = b_cds_stop
    pair.comparable_region.reverse = reverse

    logging.debug("after lcs:")
    logging.debug(pair.comparable_region)

    # TODO: fix this so glocal/global/auto are used properly
    if alignment_mode == bs_enums.ALIGNMENT_MODE.GLOBAL:
        return False

    if alignment_mode == bs_enums.ALIGNMENT_MODE.GLOCAL:
        return True

    # Region LCS: biosynthetic or min 3 domains
    if isinstance(pair.record_a, bs_gbk.Region) or isinstance(
        pair.record_b, bs_gbk.Region
    ):
        # returns True if LCS is min 3 domains, or contains a biosynthetic domain
        if len_check(pair, BigscapeConfig.REGION_MIN_LCS_LEN) or biosynthetic_check(
            pair
        ):
            return True

        # reset comparable region coordinates to full record start and stop
        logging.debug("resetting after lcs")
        reset(pair)
        return False

    else:
        # Proto LCS: must contain biosynthetic domain, no min length
        if len_check(pair, BigscapeConfig.PROTO_MIN_LCS_LEN) or biosynthetic_check(
            pair
        ):
            return True

        logging.debug("resetting after lcs")
        reset(pair)
        return False


def expand_pair(pair: RecordPair) -> bool:
    """Expand the pair

    Args:
        pair (RecordPair): pair to expand

    Returns:
        bool: True if the pair was extended, False if it does not
    """
    extend(
        pair,
        BigscapeConfig.EXPAND_MATCH_SCORE,
        BigscapeConfig.EXPAND_MISMATCH_SCORE,
        BigscapeConfig.EXPAND_GAP_SCORE,
        BigscapeConfig.EXPAND_MAX_MATCH_PERC,
    )

    # Region EXT: biosynthetic or min 5 domains
    if isinstance(pair.record_a, bs_gbk.Region) or isinstance(
        pair.record_b, bs_gbk.Region
    ):
        # returns True if LCS is min 3 domains, or contains a biosynthetic domain
        if len_check(pair, BigscapeConfig.REGION_MIN_EXPAND_LEN) or biosynthetic_check(
            pair
        ):
            return True

        # reset comparable region coordinates to full record start and stop
        logging.debug("resetting after lcs")
        reset(pair)
        return False

    # Proto EXT: min 3 domains, except no min for 1-dom rules (e.g. terpene), in which case comparable region is allowed
    # to be 1 domain as long as it is biosynthetic (already checked in do_lcs_pair)
    else:
        if (
            pair.record_a.product == pair.record_b.product
            and pair.record_a.product in BigscapeConfig.NO_MIN_CLASSES
        ):
            if len_check(pair, 0):
                return True

        else:
            if len_check(pair, BigscapeConfig.PROTO_MIN_EXPAND_LEN):
                return True

        logging.debug("resetting after extend")
        reset(pair)
        return False


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
        bs_comparison.ComparableRegion,
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
        logging.debug(pair)
        pair.log_comparable_region()
        jaccard = calc_jaccard_pair(pair, cache=False)

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
                    pair.comparable_region,
                )
            )
            continue

        # in the form [bool]. true bools means they need expansion, false they don't
        needs_expand = do_lcs_pair(pair, alignment_mode)

        if needs_expand:
            # TODO: separate these into two functions, do the extend and then calculate jaccard
            expand_pair(pair)

        if weights_label not in LEGACY_WEIGHTS:
            bin_weights = LEGACY_WEIGHTS["mix"]["weights"]
        else:
            bin_weights = LEGACY_WEIGHTS[weights_label]["weights"]
        jc_weight, ai_weight, dss_weight, anchor_boost = bin_weights

        jaccard = calc_jaccard_pair(pair)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair(pair, anchor_boost=anchor_boost)

        similarity = jaccard * jc_weight + adjacency * ai_weight + dss * dss_weight
        distance = 1 - similarity

        # at the very end, we need to inflate the comparable region coordinates to
        # include CDS without domains
        pair.comparable_region.inflate(pair)

        results.append(
            (
                pair.record_a._db_id,
                pair.record_b._db_id,
                distance,
                jaccard,
                adjacency,
                dss,
                edge_param_id,
                pair.comparable_region,
            )
        )
        logging.debug("")

    return results
