"""Contains helper functions for filtering"""

# from other modules
import big_scape.genbank as bs_gbk


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
