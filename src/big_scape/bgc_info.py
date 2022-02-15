from array import array

from src.bgctools import BgcData

class BgcInfo:
    """Class to contain info on individual BGCs"""
    name: str
    bgc_info: BgcData
    src_gbk_file: str
    
    num_genes: int
    gene_domain_counts: array
    gene_orientations: array
    bio_synth_core_positions: array

    ordered_domain_list: list
    ordered_domain_set: set

    gene_string: str
    gene_string_rev: str

    def __init__(self, name):
        self.name = name

    # dictionary of this structure:
    # BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
    #  'specific_domain_name_2'] } }
    # - cluster_name_x: cluster name (can be anything)
    # - general_domain_name_x: PFAM ID, for example 'PF00550'
    # - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names
    # in DMS unequivocally. e.g. 'PF00550_start_end', where start and end are genomic positions
    domain_name_info = {}

    def init_gene_string(self):
        """Initialize the gene string later used in DSS calculation"""
        # always find lcs seed to use for offset alignment in visualization
        # Compress the list of domains according to gene information. For example:
        # A_domlist = a b c d e f g
        # gene_domain_counts =   1  3  1  2 Number of domains per each gene in the BGC
        # gene_orientations =    1 -1 -1  1 Orientation of each gene
        # A_string = a dcb e fg List of concatenated domains
        # Takes into account gene orientation. This works effectively as putting all
        # genes in the same direction in order to be able to compare their domain content
        start = 0
        self.gene_string = []
        for num_gene in range(self.num_genes):
            domain_count = self.gene_domain_counts[num_gene]
            if self.gene_orientations[num_gene] == 1:
                # x[2:] <- small optimization, drop the "PF" from the pfam ids
                self.gene_string.append("".join(x[2:] for x in self.ordered_domain_list[start:start+domain_count]))
            else:
                a_string_doms = [self.ordered_domain_list[x][2:] for x in range(start+domain_count-1, start-1, -1)]
                self.gene_string.append("".join(a_string_doms))
            start += domain_count
        self.gene_string_rev = list(reversed(self.gene_string))
