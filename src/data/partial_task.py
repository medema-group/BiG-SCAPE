"""Contains code that scans the database for existing data and determines at which stage
of executing BiG-SCAPE the application last exited

Also contains functions to determine subsets of tasks that need to be done
"""

# from python
from __future__ import annotations
from typing import TYPE_CHECKING

# from dependencies
from sqlalchemy import select

# from other modules
from src.data import DB
import src.enums as bs_enums

# from circular imports
if TYPE_CHECKING:
    from src.genbank.gbk import GBK, CDS
    from src.hmm import HSP


def find_minimum_task(gbks: list[GBK]):
    """Finds the earliest bs_enums.TASK to start at. if new data was added, this will always
    be the load_gbks bs_enums.TASK. otherwise, it tries to find the latest bs_enums.TASK with unfinished
    business
    """
    input_data_state = get_input_data_state(gbks)

    if input_data_state.value < bs_enums.INPUT_TASK.SAME_DATA.value:
        return bs_enums.TASK.LOAD_GBKS

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

    gbk_table = DB.metadata.tables["gbk"]

    # get set of gbks in database
    db_gbk_rows = DB.execute(gbk_table.select()).all()
    db_gbk_paths = set([db_gbk_row[1] for db_gbk_row in db_gbk_rows])
    input_gbk_paths = set([str(gbk.path) for gbk in gbks])

    if db_gbk_paths == input_gbk_paths:
        return bs_enums.INPUT_TASK.SAME_DATA

    sym_dif = db_gbk_paths.symmetric_difference(input_gbk_paths)

    # still same amount in db. new data
    if len(sym_dif) == len(db_gbk_paths):
        return bs_enums.INPUT_TASK.NEW_DATA

    # same amount in new data. there was more in db than in new data
    if len(sym_dif) == len(input_gbk_paths):
        return bs_enums.INPUT_TASK.PARTIAL_DATA

    # otherwise there is some new data, some old data is missing
    return bs_enums.INPUT_TASK.MIXED_DATA


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
    scanned_cds_table = DB.metadata.tables["scanned_cds"]
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

    bgc_record_table = DB.metadata.tables["bgc_record"]

    select_statement = select(bgc_record_table.c.id)

    record_ids = set(DB.execute(select_statement).fetchall())

    distance_table = DB.metadata.tables["distance"]

    select_statement = select(distance_table.c.region_a_id).distinct()

    regions_a = set(DB.execute(select_statement).fetchall())

    select_statement = select(distance_table.c.region_b_id).distinct()

    regions_b = set(DB.execute(select_statement).fetchall())

    if len(record_ids.symmetric_difference(regions_a & regions_b)) > 0:
        return bs_enums.COMPARISON_TASK.NEW_DATA

    return bs_enums.COMPARISON_TASK.ALL_DONE
