"""Module containing helper functions to store and retrieve data used in distance calculation"""

from data.cds import gen_header
from src.big_scape.bgc_collection import BgcCollection
from src.data import get_cluster_name_list
from src.data import Database



def generate_bgc_collection(run, database: Database, BGC_INFO_DICT, GBK_FILE_DICT):
    """Generates a bgc collection for use in generate_network"""
    bgc_collection = BgcCollection()

    cluster_name_list = get_cluster_name_list(database)

    bgc_collection.initialize(cluster_name_list)
    
    rows = database.select(
        "hsp \
        join cds on hsp.cds_id = cds.id \
        join bgc on bgc.id = cds.bgc_id \
        join hmm on hmm.id = hsp.hmm_id",
        "",
        props=[
            "accession",
            "bgc.name as bgc_name"
        ]
    )

    ## prepare data
    # ordered domain list
    odomlist = dict()
    bgc_domain_name_info = dict()

    for row in rows:
        bgc_name = row["bgc_name"]
        accession = row["accession"]

        # odomlist
        if bgc_name not in odomlist:
            odomlist[bgc_name] = []
        odomlist[bgc_name].append(accession)
        
        # domain name info
        if bgc_name not in bgc_domain_name_info:
            bgc_domain_name_info[bgc_name] = {}
            if accession not in bgc_domain_name_info[bgc_name]:
                bgc_domain_name_info[bgc_name][accession] = []
            
            header = gen_header(bgc_name, row)

            bgc_domain_name_info[bgc_name][accession].append(header)
            
    ## populate the bgc collection

    # fill odomlist
    for cluster_name, domain_list in odomlist.items():
        bgc_collection.bgc_collection_dict[cluster_name].ordered_domain_list = domain_list

        # we will also want a set later on for various uses
        bgc_collection.bgc_collection_dict[cluster_name].ordered_domain_set = set(domain_list)

    bgc_collection.add_bgc_info(BGC_INFO_DICT)
    bgc_collection.add_source_gbk_files(GBK_FILE_DICT)

    # fill domain info
    for cluster_name, domain_name_info in bgc_domain_name_info.items():
        bgc_collection.bgc_collection_dict[cluster_name].domain_name_info = domain_name_info

    bgc_collection.add_gene_domain_counts(GENE_DOMAIN_COUNT)
    bgc_collection.add_bio_synth_core_pos(COREBIOSYNTHETIC_POS)
    bgc_collection.add_gene_orientations(BGC_GENE_ORIENTATION)

    # at this point we can assemble a gene string for bgc info
    bgc_collection.init_gene_strings()


    return bgc_collection

def generate_aligned_domain_seqs(run, database):
    """Generates an aligned domain seqs dictionary for use in generate_network"""
    aligned_domain_seqs = dict()
    return aligned_domain_seqs

def generate_mibig_set_indices(run, database, bgc_collection):
    """Gets the indices of mibig clusters in a bgc_collection"""
    mibig_set_indices = list()
    return mibig_set_indices

def generate_mibig_set(run, bgc_collection, mibig_set):
    """Generates a set of mibig cluster names for use in generate_network"""
    mibig_set_indices = set()
    if run.mibig.use_mibig:
        name_to_idx = {}
        for cluster_idx, cluster_name in enumerate(bgc_collection.bgc_name_tuple):
            name_to_idx[cluster_name] = cluster_idx

        for bgc in mibig_set:
            mibig_set_indices.add(name_to_idx[bgc])