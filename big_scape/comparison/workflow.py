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
import platform
from threading import Event, Condition
from typing import Generator, Callable, Optional, TypeVar, Union
from math import ceil
from .record_pair import RecordPair

# from dependencies
from sqlalchemy import select

# from other modules
from big_scape.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair
import big_scape.enums as bs_enums
import big_scape.genbank as bs_gbk
from big_scape.cli.config import BigscapeConfig
import big_scape.comparison as bs_comparison
from big_scape.data import DB
from big_scape.genbank import (
    CDS,
    GBK,
    Region,
    CandidateCluster,
    ProtoCluster,
    ProtoCore,
    BGCRecord,
)
from big_scape.hmm import HSP, HSPAlignment

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


# TODO: Test ?
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

    pair_data: Union[tuple[int, int], tuple[BGCRecord, BGCRecord]]
    if platform.system() == "Darwin":
        logging.debug(
            "Running on %s: sending full records",
            platform.system(),
        )
        pair_data = pair_generator.generate_pairs()
    else:
        logging.debug(
            "Running on %s: sending pair ids",
            platform.system(),
        )
        pair_data = pair_generator.generate_pair_ids()

    num_pairs = pair_generator.num_pairs()

    if batch_size is None:
        batch_size = get_batch_size(cores, 50000, num_pairs)
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
            # create a buffer of submitted tasks
            # to ensure that the executor is always working on something, even if
            # the main thread is busy submitting new tasks
            while len(running_futures) < max_queue_length:
                batch = batch_generator(pair_data, batch_size)

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


def do_lcs_pair(pair: RecordPair) -> bool:  # pragma no cover
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

    # Region/CandCluster LCS: biosynthetic or min 3 domains
    if (
        isinstance(pair.record_a, bs_gbk.Region)
        or isinstance(pair.record_b, bs_gbk.Region)
        or isinstance(pair.record_a, bs_gbk.CandidateCluster)
        or isinstance(pair.record_b, bs_gbk.CandidateCluster)
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

    # ProtoCluster and ProtoCore LCS: biosynthetic or min 3 domains
    else:
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

    # Region/CandCluster EXT checks: biosynthetic and min 3 or min 5 domains,
    # except no min for 1-dom rules (e.g. terpene)
    if (
        isinstance(pair.record_a, bs_gbk.Region)
        or isinstance(pair.record_b, bs_gbk.Region)
        or isinstance(pair.record_a, bs_gbk.CandidateCluster)
        or isinstance(pair.record_b, bs_gbk.CandidateCluster)
    ):
        # returns True if EXT is min 5 domains, min 0 domains and contains biosynthetic domain for no_min_class
        if (
            pair.record_a.product == pair.record_b.product
            and pair.record_a.product in BigscapeConfig.NO_MIN_CLASSES
        ):
            if len_check(
                pair, BigscapeConfig.REGION_MIN_EXPAND_LEN
            ) or biosynthetic_check(pair):
                return True

        # returns True if EXT is min 5 domains, or contains a biosynthetic domain and is min 3 domains
        else:
            if len_check(pair, BigscapeConfig.REGION_MIN_EXPAND_LEN) or (
                len_check(pair, BigscapeConfig.REGION_MIN_EXPAND_LEN_BIO)
                and biosynthetic_check(pair)
            ):
                return True

    # ProtoCluster EXT: min 3 domains and biosynthetic, except no min for 1-dom rules (e.g. terpene)
    elif isinstance(pair.record_a, bs_gbk.ProtoCluster) and isinstance(
        pair.record_b, bs_gbk.ProtoCluster
    ):
        # no_min_class: biosynthetic & min_len = 0
        if (
            pair.record_a.product == pair.record_b.product
            and pair.record_a.product in BigscapeConfig.NO_MIN_CLASSES
        ):
            if biosynthetic_check(pair):
                return True

        # biosynthetic & min_len = 3
        else:
            if len_check(
                pair, BigscapeConfig.PROTO_MIN_EXPAND_LEN
            ) and biosynthetic_check(pair):
                return True

    # ProtoCore Ext: needs to be biosynthetic
    elif isinstance(pair.record_a, bs_gbk.ProtoCore) and isinstance(
        pair.record_b, bs_gbk.ProtoCore
    ):
        if biosynthetic_check(pair):
            return True

    else:
        logging.debug(
            "comparing protocluster to protocore, something must have gone wrong here, pair: %s",
            pair,
        )
        reset(pair)
        return False

    logging.debug("resetting after extend")
    reset(pair)
    return False


# TODO: Test ?
def calculate_scores_pair(
    data: tuple[
        list[Union[tuple[int, int], tuple[BGCRecord, BGCRecord]]],
        bs_enums.ALIGNMENT_MODE,
        int,
        str,
    ]
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
        data (tuple[list[tuple[int, int]], str, str]): list of pairs, alignment mode,
        bin label

    Returns:
        list[tuple[int, int, float, float, float, float, int, int, int, int, int, int,
        int, int, bool, str,]]: list of scores for each pair in the
        order as the input data list, including lcs and extension coordinates
    """
    data, alignment_mode, edge_param_id, weights_label = data

    # convert database ids to minimal record objects
    if isinstance(data[0][0], int):
        pair_ids = data
        records = fetch_records_from_database(pair_ids)
    else:
        for pair in data:
            records = {record._db_id: record for record in pair}
            pair_ids = [
                (record.record_a._db_id, record.record_b._db_id) for record in pair
            ]

    results = []

    # TODO: this fails since DB getting accessed from child processes
    # seems to be a problem with the DB connection (for mac?)
    # weights_label = bs_comparison.get_edge_weight(edge_param_id)

    for id_a, id_b in pair_ids:
        if id_a not in records or id_b not in records:
            comparable_region = bs_comparison.ComparableRegion(
                0, 0, 0, 0, 0, 0, 0, 0, False
            )
            results.append(
                (id_a, id_b, 1.0, 0.0, 0.0, 0.0, edge_param_id, comparable_region)
            )
            continue

        pair = RecordPair(records[id_a], records[id_b])
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

        # GLOBAL/GLOCAL/AUTO
        # GLOBAL the comparable region is the full record
        # GLOCAL computes a the comparable region based on LCS EXT
        # AUTO computes a the comparable region based on LCS EXT if
        # either record is on a contig edge

        if alignment_mode == bs_enums.ALIGNMENT_MODE.GLOCAL or (
            alignment_mode == bs_enums.ALIGNMENT_MODE.AUTO
            and (pair.record_a.contig_edge or pair.record_b.contig_edge)
        ):
            needs_expand = do_lcs_pair(pair)
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


def fetch_records_from_database(pairs: list[tuple[int, int]]) -> dict[int, BGCRecord]:
    """Constructs an index containing minimal BGCRecord objects. This can be used from
    threads to gather the minimal information needed to perform distance calculations.

    Args:
        pairs (list[tuple[int, int]]): list of record pair in database ids

    Returns:
        dict[int, BGCRecord]: index linking each database id to a minimal record object
    """

    # First gather records to collect from the database
    pair_ids = set([db_id for pair in pairs for db_id in pair])

    if DB.metadata is None:
        raise RuntimeError("DB metadata is None!")

    gbk_table = DB.metadata.tables["gbk"]
    record_table = DB.metadata.tables["bgc_record"]
    cds_table = DB.metadata.tables["cds"]
    hsp_table = DB.metadata.tables["hsp"]
    algn_table = DB.metadata.tables["hsp_alignment"]

    # gather minimally needed information for distance calculation and object creation
    query_statement = (
        select(
            gbk_table.c.id,
            record_table.c.id,
            record_table.c.record_type,
            record_table.c.nt_start,
            record_table.c.nt_stop,
            record_table.c.contig_edge,
            cds_table.c.id,
            cds_table.c.nt_start,
            cds_table.c.nt_stop,
            cds_table.c.orf_num,
            cds_table.c.strand,
            cds_table.c.gene_kind,
            hsp_table.c.id,
            hsp_table.c.accession,
            hsp_table.c.bit_score,
            hsp_table.c.env_start,
            hsp_table.c.env_stop,
            algn_table.c.alignment,
        )
        .join(record_table, record_table.c.gbk_id == gbk_table.c.id)
        .join(cds_table, cds_table.c.gbk_id == gbk_table.c.id)
        .join(hsp_table, hsp_table.c.cds_id == cds_table.c.id)
        .join(algn_table, algn_table.c.hsp_id == hsp_table.c.id)
        .where(record_table.c.id.in_(pair_ids))
    )

    pair_data = DB.execute(query_statement).fetchall()

    if not pair_data:
        raise RuntimeError("Data of pairs not found in database!")

    # make an index of formatted data for each record to be reused
    gbk_index: dict[int, GBK] = {}
    record_index: dict[int, BGCRecord] = {}
    cds_index: dict[int, CDS] = {}
    added_hsp_ids: set[int] = set()

    # iteratively build the required data structures
    for row in pair_data:
        gbk_id, rec_id, cds_id, hsp_id = row[0], row[1], row[6], row[12]

        # create dummy GBK object to link records and CDSs
        if gbk_id not in gbk_index:
            gbk_index[gbk_id] = GBK("", "", "")

        # keep track of record_type to enable comparison on e.g. protocluster level
        if rec_id not in record_index:
            record_type = row[2]
            if record_type == "region":
                record_index[rec_id] = Region(
                    gbk_index[gbk_id], 0, row[3], row[4], row[5], ""
                )
            elif record_type == "cand_cluster":
                record_index[rec_id] = CandidateCluster(
                    gbk_index[gbk_id], 0, row[3], row[4], row[5], "", "", {}
                )
            elif record_type == "protocluster":
                record_index[rec_id] = ProtoCluster(
                    gbk_index[gbk_id], 0, row[3], row[4], row[5], "", {}
                )
            elif record_type == "proto_core":
                record_index[rec_id] = ProtoCore(
                    gbk_index[gbk_id], 0, row[3], row[4], row[5], ""
                )
            record_index[rec_id]._db_id = rec_id

        # simultaneously save CDSs and assign to the correct parent GBK
        if cds_id not in cds_index:
            cds_index[cds_id] = CDS(row[7], row[8])
            cds_index[cds_id].orf_num = row[9]
            cds_index[cds_id].strand = row[10]
            cds_index[cds_id].gene_kind = row[11]
            gbk_index[gbk_id].genes.append(cds_index[cds_id])

        # finally, assign HSPs to their respective CDS
        if hsp_id not in added_hsp_ids:
            new_hsp = HSP(cds_index[cds_id], row[13], row[14], row[15], row[16])
            new_hsp.alignment = HSPAlignment(new_hsp, row[17])
            cds_index[cds_id].hsps.append(new_hsp)
            added_hsp_ids.add(hsp_id)

    return record_index
