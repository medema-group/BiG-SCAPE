"""Module containing a class to store BGC info when read from the legacy GBK
reader

Author: Jorge Navarro, Arjan Draisma"""


class BgcData:
    """Class to contain BGC data"""
    def __init__(self, accession_id, description, product, records, max_width, bgc_size, organism,
                 taxonomy, biosynthetic_genes, contig_edge):
        # These two properties come from the genbank file:
        self.accession_id = accession_id
        self.description = description
        # AntiSMASH predicted class of compound:
        self.product = product
        # number of records in the genbank file (think of multi-locus BGCs):
        self.records = records
        # length of largest record (it will be used for ArrowerSVG):
        self.max_width = int(max_width)
        # length of the entire bgc (can include several records/subclusters)
        self.bgc_size = bgc_size
        # organism
        self.organism = organism
        # taxonomy as a string (of comma-separated values)
        self.taxonomy = taxonomy
        # Internal set of tags corresponding to genes that AntiSMASH marked
        # as "Kind: Biosynthetic". It is formed as
        # clusterName + "_ORF" + cds_number + ":gid:" + gene_id + ":pid:" + protein_id + ":loc:"
        # + gene_start + ":" + gene_end + ":strand:" + {+,-}
        self.biosynthetic_genes = biosynthetic_genes
        # AntiSMASH 4+ marks BGCs that sit on the edge of a contig
        self.contig_edge = contig_edge
