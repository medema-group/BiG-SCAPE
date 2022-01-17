from collections import defaultdict
from xml import dom


class ClusterInfo:
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

    dom_start: int
    dom_end: int

    gene_string: str

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

        self.gene_string = []

    def init_dom_seq_slices(self, bgcs):
        # initialize domain sequence slices
        # They might change if we manage to find a valid overlap
        for domain in self.domlist_set:
            self.domain_seq_slice_bottom[domain] = 0
            self.domain_seq_slice_top[domain] = len(bgcs[self.cluster_name][domain])
    
    def init_dom_borders(self):
        self.dom_start = 0
        self.dom_end = len(self.domlist)
    
    def init_gene_string(self):
        # Compress the list of domains according to gene information. For example:
        # A_domlist = a b c d e f g
        # cluster_info_a.dcg =   1  3  1  2 Number of domains per each gene in the BGC
        # cluster_info_a.go =    1 -1 -1  1 Orientation of each gene
        # A_string = a dcb e fg List of concatenated domains
        # Takes into account gene orientation. This works effectively as putting all
        # genes in the same direction in order to be able to compare their domain content
        start = 0
        for num_gene in range(self.num_genes):
            domain_count = self.dcg[num_gene]
            if self.go[num_gene] == 1:
                # x[2:] <- small optimization, drop the "PF" from the pfam ids
                self.gene_string.append("".join(x[2:] for x in self.domlist[start:start+domain_count]))
            else:
                a_string_doms = [self.domlist[x][2:] for x in range(start+domain_count-1, start-1, -1)]
                self.gene_string.append("".join(a_string_doms))
            start += domain_count

