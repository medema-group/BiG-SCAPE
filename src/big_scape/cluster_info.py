from collections import defaultdict
from xml import dom


class cluster_info:
    # base info
    cluster_name: str
    domlist: list
    dcg: list
    core_pos: any
    go: any

    # calculated afterwards
    # Number of genes in each BGC
    num_genes: int

    domlist_set: set

    # define the subset of domain sequence tags to include in
    # the DSS calculation. This is done for every domain.
    domain_seq_slice_top: dict
    domain_seq_slice_bottom: dict

    num_anchor_domains: int
    num_non_anchor_domains: int

    def __init__(self, cluster_name, domlist, dcg, core_pos, go):
        self.cluster_name = cluster_name
        self.domlist = domlist
        self.dcg = dcg
        self.core_pos = core_pos
        self.go = go
        self.initialize()

    def initialize(self):
        self.num_genes = len(self.dcg)
        self.domlist_set = set(self.domlist)

        self.domain_seq_slice_bottom = defaultdict(int)
        self.domain_seq_slice_top = defaultdict(int)

    def init_dom_seq_slices(self, bgcs):
        # initialize domain sequence slices
        # They might change if we manage to find a valid overlap
        for domain in self.domlist_set:
            self.domain_seq_slice_bottom[domain] = 0
            self.domain_seq_slice_top[domain] = len(bgcs[self.cluster_name][domain])

