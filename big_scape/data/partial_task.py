"""Contains code that scans the database for existing data and determines at which stage
of executing BiG-SCAPE the application last exited

Also contains functions to determine subsets of tasks that need to be done
"""

# from python
from __future__ import annotations
from typing import TYPE_CHECKING, Generator, Optional

# from dependencies
from sqlalchemy import select

# from other modules
from big_scape.data import DB
import big_scape.enums as bs_enums

# from circular imports
if TYPE_CHECKING:
    from big_scape.comparison import RecordPairGenerator
    from big_scape.genbank.gbk import GBK, CDS
    from big_scape.hmm import HSP


def find_minimum_task(gbks: list[GBK]):
    """Finds the earliest bs_enums.TASK to start at. if new data was added, this will
    always be the load_gbks bs_enums.TASK. otherwise, it tries to find the latest
    bs_enums.TASK with unfinished business
    """
    input_data_state = get_input_data_state(gbks)

    # new data or mixed data
    if (
        input_data_state.value == bs_enums.INPUT_TASK.NEW_DATA.value
        or input_data_state.value == bs_enums.INPUT_TASK.MIXED_DATA.value
        or input_data_state.value == bs_enums.INPUT_TASK.NO_DATA.value
    ):
        # gbks from input need to be loaded into the in-memory database
        return bs_enums.TASK.SAVE_GBKS

    hmm_data_state = get_hmm_data_state(gbks)

    if hmm_data_state.value < bs_enums.HMM_TASK.NEED_ALIGN.value:
        return bs_enums.TASK.HMM_SCAN

    if hmm_data_state.value < bs_enums.HMM_TASK.ALL_ALIGNED.value:
        return bs_enums.TASK.HMM_ALIGN

    comparison_data_state = get_comparison_data_state(gbks)

    if comparison_data_state.value < bs_enums.COMPARISON_TASK.ALL_DONE.value:
        return bs_enums.TASK.COMPARISON

    return bs_enums.TASK.NOTHING_TO_DO


def get_input_data_state(gbks: list[GBK]) -> bs_enums.INPUT_TASK:
    """Returns the status of input data (gbks and regions) in the in-memory database"""
    distance_count = DB.get_table_row_count("gbk")

    if distance_count == 0:
        return bs_enums.INPUT_TASK.NO_DATA

    gbk_table = DB.get_table("gbk")

    # get set of gbks in database
    db_gbk_rows = DB.execute(gbk_table.select()).all()
    db_gbk_hashes: set[str] = {db_gbk_row[2] for db_gbk_row in db_gbk_rows}
    input_gbk_hashes: set[str] = {str(gbk.hash) for gbk in gbks}

    if db_gbk_hashes == input_gbk_hashes:
        return bs_enums.INPUT_TASK.SAME_DATA

    union = db_gbk_hashes & input_gbk_hashes

    # all new data
    if len(union) == 0:
        return bs_enums.INPUT_TASK.NEW_DATA

    # only partial data which is already in database
    if len(union) == len(input_gbk_hashes):
        return bs_enums.INPUT_TASK.PARTIAL_DATA

    # otherwise there is some new data, some old data is missing
    return bs_enums.INPUT_TASK.MIXED_DATA


def get_missing_gbks(gbks: list[GBK]) -> list[GBK]:
    """Find which GBKs are missing from the database and return them

    Args:
        gbks (list[GBK]): List of GBKs to check

    Returns:
        list[GBK]: List of GBKs that are missing from the database
    """
    # dictionary of gbk path to gbk object
    gbk_dict = {str(gbk.hash): gbk for gbk in gbks}

    gbk_table = DB.get_table("gbk")

    # get set of gbks in database
    db_gbk_rows = DB.execute(gbk_table.select()).all()
    db_gbk_hashes: set[int] = {db_gbk_row[2] for db_gbk_row in db_gbk_rows}

    missing_gbks = []

    for gbk_hash in gbk_dict:
        if gbk_hash not in db_gbk_hashes:
            missing_gbks.append(gbk_dict[gbk_hash])

    return missing_gbks


def get_hmm_data_state(gbks: list[GBK]) -> bs_enums.HMM_TASK:
    """Retuns the state of data hmm processing in the in-memory database"""
    hsp_count = DB.get_table_row_count("hsp")

    if hsp_count == 0:
        return bs_enums.HMM_TASK.NO_DATA

    cds_count = DB.get_table_row_count("cds")
    scanned_cds_count = DB.get_table_row_count("scanned_cds")

    if cds_count > scanned_cds_count:
        return bs_enums.HMM_TASK.NEED_SCAN

    align_count = DB.get_table_row_count("hsp_alignment")

    if align_count < hsp_count:
        return bs_enums.HMM_TASK.NEED_ALIGN

    return bs_enums.HMM_TASK.ALL_ALIGNED


def get_cds_to_scan(gbks: list[GBK]) -> list[CDS]:
    """Find which cds within a list of GBK objects need to be scanned using hmmscan

    Args:
        gbks (list[GBK]): List of GBKs

    Returns:
        list[CDS]: List of CDS which were not scanned according to the data in scanned_cds table
    """
    hmm_state = get_hmm_data_state(gbks)

    if hmm_state == bs_enums.HMM_TASK.NO_DATA:
        cds_to_scan = []
        for gbk in gbks:
            cds_to_scan.extend(gbk.genes)
        return cds_to_scan

    cds_to_scan = []

    # get a list of database cds_ids that are present in the cds_scanned table

    scanned_cds_table = DB.get_table("scanned_cds")
    select_query = select(scanned_cds_table.c.cds_id)
    scanned_cds_ids = set(DB.execute(select_query))

    for gbk in gbks:
        for gene in gbk.genes:
            if (gene._db_id,) not in scanned_cds_ids:
                cds_to_scan.append(gene)

    return cds_to_scan


def get_hsp_to_align(gbks: list[GBK]) -> list[HSP]:
    """Get a list of HSP objects that lack an alignment

    Args:
        gbks (list[GBK]): list of GBKs tos earch for missing gene hsp aligments

    Returns:
        list[HSP]: List of HSPs with missing alignments
    """
    hsps_to_align = []
    for gbk in gbks:
        for gene in gbk.genes:
            for hsp in gene.hsps:
                if hsp.alignment is None:
                    hsps_to_align.append(hsp)
    return hsps_to_align


def get_comparison_data_state(gbks: list[GBK]) -> bs_enums.COMPARISON_TASK:
    """Retuns the state of pairwise comparison data in the in-memory database"""

    distance_count = DB.get_table_row_count("distance")

    if distance_count == 0:
        return bs_enums.COMPARISON_TASK.NO_DATA

    # TODO: this needs changing once we implement protocluster/protocore
    # stuff. currently this is a naive way of calculating this

    # check if all record ids are present in the comparison region ids

    bgc_record_table = DB.get_table("bgc_record")

    select_statement = select(bgc_record_table.c.id)

    record_ids = set(DB.execute(select_statement).fetchall())

    distance_table = DB.get_table("distance")

    select_statement = select(distance_table.c.record_a_id).distinct()

    regions_a = set(DB.execute(select_statement).fetchall())

    select_statement = select(distance_table.c.record_b_id).distinct()

    regions_b = set(DB.execute(select_statement).fetchall())

    if len(record_ids.symmetric_difference(regions_a & regions_b)) > 0:
        return bs_enums.COMPARISON_TASK.NEW_DATA

    return bs_enums.COMPARISON_TASK.ALL_DONE


# TODO: does not seem to be used
def get_missing_distances(
    pair_generator: RecordPairGenerator,
) -> Generator[tuple[Optional[int], Optional[int]], None, None]:
    """Get a generator of BGCPairs that are missing from a network

    Args:
        network (BSNetwork): network to check
        bin (BGCBin): bin to check

    Yields:
        Generator[BGCPair]: generator of BGCPairs that are missing from the network
    """

    distance_table = DB.get_table("distance")

    # get all region._db_id in the bin
    select_statement = (
        select(distance_table.c.record_a_id, distance_table.c.record_b_id)
        .where(distance_table.c.record_a_id.in_(pair_generator.record_ids))
        .where(distance_table.c.record_b_id.in_(pair_generator.record_ids))
    )

    # generate a set of tuples of region id pairs
    existing_distances = set(DB.execute(select_statement).fetchall())

    for pair in pair_generator.generate_pair_ids():
        # if the pair is not in the set of existing distances, yield it
        if pair not in existing_distances and pair[::-1] not in existing_distances:
            yield pair
