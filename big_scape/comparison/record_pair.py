from big_scape.comparison.comparable_region import ComparableRegion
from big_scape.genbank import BGCRecord
from big_scape.hmm import HSP


import logging
from typing import Optional


class RecordPair:
    """Contains a pair of BGC records, which can be any type of BGCRecord

    This will also contain any other necessary information specific to this pair needed
    to generate the scores

    Attributes:
        record_a (BGCRecord): First BGC record in the pair
        record_b (BGCRecord): Second BGC record in the pair
        comparable_region (ComparableRegion): A comparable region between the two records

        domain_lists: Optional[tuple[list[HSP], list[HSP]]]
        domain_sets: Optional[tuple[set[HSP], set[HSP]]]
        domain_dicts: Optional[tuple[dict[HSP, list[int]], dict[HSP, list[int]]]]
    """

    def __init__(self, record_a: BGCRecord, record_b: BGCRecord):
        self.record_a = record_a
        self.record_b = record_b

        if record_a.parent_gbk is None or record_b.parent_gbk is None:
            raise ValueError("Region in pair has no parent GBK!")

        # comparable regions start "deflated", meaning only CDS with domains
        a_len = len(record_a.get_cds_with_domains())
        b_len = len(record_b.get_cds_with_domains())
        a_domain_len = len(record_a.get_hsps())
        b_domain_len = len(record_b.get_hsps())

        self.comparable_region: ComparableRegion = ComparableRegion(
            0, a_len, 0, b_len, 0, a_domain_len, 0, b_domain_len, False
        )

        self.domain_lists: Optional[tuple[list[HSP], list[HSP]]] = None
        self.domain_sets: Optional[tuple[set[HSP], set[HSP]]] = None
        self.domain_dicts: Optional[
            tuple[dict[HSP, list[int]], dict[HSP, list[int]]]
        ] = None

    def __repr__(self) -> str:
        return f"Pair {self.record_a} - {self.record_b}"

    def __hash__(self) -> int:
        a_hash = hash(self.record_a)
        b_hash = hash(self.record_b)

        # order doesn't matter
        return a_hash + b_hash

    def __eq__(self, _o) -> bool:
        if not isinstance(_o, RecordPair):
            return False

        if self.record_a == _o.record_a and self.record_b == _o.record_b:
            return True
        if self.record_a == _o.record_b and self.record_b == _o.record_a:
            return True

        return False

    def get_domain_sets(
        self, regenerate=False, cache=True
    ) -> tuple[set[HSP], set[HSP]]:
        """Returns a tuple containing sets of domains within the comparable region of
        two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True. This is done to minimize the computational cost of
        intializing a set from a list at the cost of memory.

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one
            cache (bool): whether to cache the result of this operation

        Returns:
            tuple[set[HSP], set[HSP]]
        """

        if regenerate or self.domain_sets is None:
            a_domain_list, b_domain_list = self.get_domain_lists(cache=cache)

            if not cache:
                return (set(a_domain_list), set(b_domain_list))

            self.domain_sets = (set(a_domain_list), set(b_domain_list))

        return self.domain_sets

    def get_domain_lists(
        self, regenerate=False, cache=True, reverse=None
    ) -> tuple[list[HSP], list[HSP]]:
        """Returns a tuple corresponding (ordered) lists of CDS domains within the
        comparable region of two BGCs

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True. This is done to minimize the computational cost of
        intializing a set from a list at the cost of memory.

        If reverse is set to none, uses the reverse propoerty on comparable region

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one,
            cache (bool): whether to cache the result for faster retrieval
            reverse (bool, optional): Whether to return a reversed list for region B.
            Defaults to None.

        Returns:
            tuple[list[HSP], list[HSP]]
        """

        if regenerate or self.domain_lists is None:
            a_start = self.comparable_region.a_start
            a_stop = self.comparable_region.a_stop

            b_start = self.comparable_region.b_start
            b_stop = self.comparable_region.b_stop

            reverse = self.comparable_region.reverse

            a_cds_list = self.record_a.get_cds_with_domains()[a_start:a_stop]
            b_cds_list = self.record_b.get_cds_with_domains(reverse=reverse)[
                b_start:b_stop
            ]

            a_domain_list: list[HSP] = []
            for a_cds in a_cds_list:
                a_domain_list.extend(a_cds.hsps)

            b_domain_list: list[HSP] = []
            for b_cds in b_cds_list:
                b_domain_list.extend(b_cds.hsps)

            if not cache:
                return (a_domain_list, b_domain_list)

            self.domain_lists = (a_domain_list, b_domain_list)

        return self.domain_lists

    def get_domain_dicts(
        self, regenerate=False
    ) -> tuple[dict[HSP, list[int]], dict[HSP, list[int]]]:
        """Returns a dictionary of domains for each BGC in this comaprable region.
        Dictionary keys are domains and ints are the index of that domain in the domain
        list from which the dictionary was generated. This will always be the list
        returned by self.get_domain_list at the moment this method is called

        This method caches the result and re-uses the result on further calls unless
        regenerate is set to True

        Args:
            regenerate (bool): whether to replace the cached result set with a newly
            generated one

        Returns:
            tuple[dict[HSP, list[int]], dict[HSP, list[int]]]: dictionary of domain
            accessions to list indexes
        """

        if regenerate or self.domain_dicts is None:
            domain_dict_a: dict[HSP, list[int]] = {}
            domain_dict_b: dict[HSP, list[int]] = {}

            domain_list_a, domain_list_b = self.get_domain_lists()

            for idx_a, domain_a in enumerate(domain_list_a):
                if domain_a not in domain_dict_a:
                    domain_dict_a[domain_a] = []

                domain_dict_a[domain_a].append(idx_a)

            for idx_b, domain_b in enumerate(domain_list_b):
                if domain_b not in domain_dict_b:
                    domain_dict_b[domain_b] = []

                domain_dict_b[domain_b].append(idx_b)

            self.domain_dicts = (domain_dict_a, domain_dict_b)

        return self.domain_dicts

    def log_comparable_region(self, label="<") -> None:  # pragma: no cover
        """Prints a debug level log of the comparable region

        Args:
            label (str): A string to visually indicate the comparable region
        """
        if logging.getLogger().level > logging.DEBUG:
            return

        a_cds_list = self.record_a.get_cds_with_domains()
        b_cds_list = self.record_b.get_cds_with_domains(
            reverse=self.comparable_region.reverse
        )

        b_start = self.comparable_region.b_start
        b_stop = self.comparable_region.b_stop

        for i in range(max(len(a_cds_list), len(b_cds_list))):
            a_domains = ""
            if i < len(a_cds_list):
                a_domains = "|".join(
                    [hsp.domain[2:].split(".")[0] for hsp in a_cds_list[i].hsps]
                )

            b_domains = ""
            if i < len(b_cds_list):
                b_domains = "|".join(
                    [hsp.domain[2:].split(".")[0] for hsp in b_cds_list[i].hsps]
                )

            a_region_str = ""
            b_region_str = ""

            a_in_region = (
                i > self.comparable_region.a_start and i < self.comparable_region.a_stop
            )
            if a_in_region:
                a_region_str = label
            if i == self.comparable_region.a_start:
                a_region_str = "START"
            if i == self.comparable_region.a_stop:
                a_region_str = "STOP"

            b_in_region = i > b_start and i < b_stop
            if b_in_region:
                b_region_str = label
            if i == b_start:
                b_region_str = "START"
            if i == b_stop:
                b_region_str = "STOP"

            log_line = " ".join(
                [
                    f"{i:<3}",
                    f"{a_domains:<25}",
                    f"{a_region_str:<10}",
                    f"{b_domains:<25}",
                    f"{b_region_str:<10}",
                ]
            )
            logging.debug(log_line)
