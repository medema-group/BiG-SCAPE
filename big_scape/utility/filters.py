"""Contains helper functions for filtering"""

import logging

# from other modules
import big_scape.genbank as bs_gbk
from big_scape.comparison import get_record_category


def domain_includelist_filter(run: dict, all_bgc_records: list[bs_gbk.BGCRecord]):
    include_domain_accessions = []
    domainlist_bgc_records = []

    if run["domain_includelist_any"] is not None:
        include_domain_accessions = run["domain_includelist_any"]

        for record in all_bgc_records:
            record_hsps = record.get_hsps()

            result = any(hsp.domain in include_domain_accessions for hsp in record_hsps)

            if result:
                domainlist_bgc_records.append(record)

        return domainlist_bgc_records

    if run["domain_includelist_all"] is not None:
        include_domain_accessions = run["domain_includelist_all"]

        for record in all_bgc_records:
            record_domains = set(hsp.domain for hsp in record.get_hsps())
            # TODO: is this more or less efficient than using a record_domains generator?

            result = all(
                domain in record_domains for domain in include_domain_accessions
            )

            if result:
                domainlist_bgc_records.append(record)

        return domainlist_bgc_records


def generic_category_class_filter(
    class_categ: set[str], include: set[str], exclude: set[str]
) -> bool:
    """Generic function for both class and category based filtering of a single record

    exclude class/category is checked first, failing any records that contain any of
    these classes/categories. Then, include class/category is checked, passing only any
    records that contain any of these classes/categories.

    Args:
        class_categ (set[str]): set of classes or categories present in a single record
        include (set[str]): set of classes or categories to include
        exclude (set[str]): set of classes or categories to exclude

    Returns:
        bool: True if the current record passes filters, False if it is filtered out
    """
    if exclude:
        if class_categ & exclude:
            return False

    if include:
        if class_categ & include:
            return True
        return False

    return True


def category_filter(
    run: dict, all_bgc_records: list[bs_gbk.BGCRecord]
) -> list[bs_gbk.BGCRecord]:
    """Filter working records based on categories

    Args:
        run (dict): run parameters
        all_bgc_records (list[bs_gbk.BGCRecord]): working records

    Raises:
        RuntimeError: All records were filtered out

    Returns:
        list[bs_gbk.BGCRecord]: list of records passing set filters
    """
    orig_size = len(all_bgc_records)
    include_categs = run["include_categories"]
    exclude_categs = run["exclude_categories"]

    filtered_bgc_records = []
    for record in all_bgc_records:
        record_categs = set(get_record_category(record).split("."))
        if generic_category_class_filter(record_categs, include_categs, exclude_categs):
            filtered_bgc_records.append(record)

    if len(filtered_bgc_records) == 0:
        logging.error(
            "No BGC records remain after include/exclude categories filtering",
        )
        raise RuntimeError()

    if orig_size != len(filtered_bgc_records):
        logging.info(
            "Continuing with %s BGC records after include/exclude categories filtering",
            len(filtered_bgc_records),
        )
    return filtered_bgc_records


def class_filter(
    run: dict, all_bgc_records: list[bs_gbk.BGCRecord]
) -> list[bs_gbk.BGCRecord]:
    """Filter working records based on classes

    Args:
        run (dict): run parameters
        all_bgc_records (list[bs_gbk.BGCRecord]): working records

    Raises:
        RuntimeError: All records were filtered out

    Returns:
        list[bs_gbk.BGCRecord]: list of records passing set filters
    """
    orig_size = len(all_bgc_records)
    include_class = run["include_classes"]
    exclude_class = run["exclude_classes"]

    filtered_bgc_records = []
    for record in all_bgc_records:
        rec_classes = set(record.product.split("."))
        if generic_category_class_filter(rec_classes, include_class, exclude_class):
            filtered_bgc_records.append(record)

    if len(filtered_bgc_records) == 0:
        logging.error(
            "No BGC records remain after include/exclude classes filtering",
        )
        raise RuntimeError()

    if orig_size != len(filtered_bgc_records):
        logging.info(
            "Continuing with %s BGC records after include/exclude classes filtering",
            len(filtered_bgc_records),
        )
    return filtered_bgc_records
