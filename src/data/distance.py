"""Module containing helper functions to store and retrieve data used in distance calculation"""

from array import array
from src.data.cds import gen_header
from src.big_scape.bgc_collection import BgcCollection
from src.data import get_cluster_name_list
from src.data import Database



def generate_bgc_collection(run, database: Database, BGC_INFO_DICT, GBK_FILE_DICT):
    """Generates a bgc collection for use in generate_network"""
    bgc_collection = BgcCollection()

    cluster_name_list = get_cluster_name_list(database)

    bgc_collection.initialize(cluster_name_list)

    # protein domain needs to be sorted by the absolute start
    # we could do this in python but somehow it seems it might be easier to do in
    # sqlite
    # the CASE statement here calculates the absolute start of the domain
    # this is used to order by absolute start
    rows = database.select(
        "hsp \
        JOIN hsp_alignment ON hsp_alignment.hsp_id = hsp.id \
        JOIN cds ON hsp.cds_id = cds.id \
        JOIN bgc ON bgc.id = cds.bgc_id \
        JOIN hmm ON hmm.id = hsp.hmm_id",
        "ORDER BY bgc_id, absolute_start",
        props=[
            "accession",
            "bgc.name as bgc_name",
            "nt_start",
            "nt_end",
            "strand",
            "env_start",
            "env_end",
            "orf_id",
            "CASE WHEN strand = 1 \
                THEN 3 * env_start + nt_start \
            ELSE \
                nt_end - 3 * env_start - 3 * (env_end - env_start) \
            END AS absolute_start"
        ]
    )

    ## prepare data
    # ordered domain list and bgc domain name info
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

        # assemble header expected later on
        header = gen_header(bgc_name, row)[1:] + ":" + str(row["env_start"]) + ":" + str(row["env_end"])

        bgc_domain_name_info[bgc_name][accession].append(header)

    # domain count, core biosynth pos and gene orientation
    gene_domain_count = dict()
    corebiosynthetic_pos = dict()
    bgc_gene_orientation = dict()

    rows = database.select(
        "hsp \
        join cds on hsp.cds_id = cds.id \
        join bgc on bgc.id = cds.bgc_id \
        join hmm on hmm.id = hsp.hmm_id",
        "group by cds_id",
        props=[
            "bgc.name as bgc_name",
            "count(cds_id) as cds_count",
            "biosynthetic",
            "strand",
            "orf_id"
        ]
    )


    for idx, row in enumerate(rows):
        bgc_name = row["bgc_name"]
        if bgc_name not in gene_domain_count:
            gene_domain_count[bgc_name] = array('B')
        if bgc_name not in corebiosynthetic_pos:
            corebiosynthetic_pos[bgc_name] = array('H')
        if bgc_name not in bgc_gene_orientation:
            bgc_gene_orientation[bgc_name] = array('b')

        gene_domain_count[bgc_name].append(row["cds_count"])
        if row["biosynthetic"] == 1:
            corebiosynthetic_pos[bgc_name].append(row["orf_id"])
        bgc_gene_orientation[bgc_name].append(row["strand"])

    ## populate the bgc collection

    # fill odomlist
    bgc_collection.bgc_ordered_domain_list = odomlist
    for cluster_name, domain_list in odomlist.items():
        bgc_collection.bgc_collection_dict[cluster_name].ordered_domain_list = domain_list

        # we will also want a set later on for various uses
        bgc_collection.bgc_collection_dict[cluster_name].ordered_domain_set = set(domain_list)

    bgc_collection.add_bgc_info(BGC_INFO_DICT)
    bgc_collection.add_source_gbk_files(GBK_FILE_DICT)

    # fill domain info
    for cluster_name, domain_name_info in bgc_domain_name_info.items():
        bgc_collection.bgc_collection_dict[cluster_name].domain_name_info = domain_name_info

    bgc_collection.add_gene_domain_counts(gene_domain_count)
    bgc_collection.add_bio_synth_core_pos(corebiosynthetic_pos)
    bgc_collection.add_gene_orientations(bgc_gene_orientation)

    # at this point we can assemble a gene string for bgc info
    bgc_collection.init_gene_strings()


    return bgc_collection

def generate_aligned_domain_seqs(run, database):
    """Generates an aligned domain seqs dictionary for use in generate_network"""
    aligned_domain_seqs = dict()

    rows = database.select(
        "hsp \
        join hsp_alignment on hsp_alignment.hsp_id = hsp.id \
        join cds on hsp.cds_id = cds.id \
        join bgc on bgc.id = cds.bgc_id \
        join msa on msa.cds_id = hsp.cds_id and msa.hmm_id = hsp.hmm_id",
        "",
        props=[
            "bgc.name as bgc_name",
            "nt_start",
            "nt_end",
            "msa.env_start",
            "msa.env_end",
            "strand",
            "algn_string",
            "orf_id"
        ]
    )

    for row in rows:
        bgc_name = row["bgc_name"]
        header = gen_header(bgc_name, row)[1:] + ":" + str(row["env_start"]) + ":" + str(row["env_end"])
        aligned_domain_seqs[header] = row["algn_string"]


    return aligned_domain_seqs

def generate_mibig_set_indices(run, bgc_collection: BgcCollection, mibig_set):
    """Gets the indices of mibig clusters in a bgc_collection and mibig_set"""
    mibig_set_indices = set()
    if run.mibig.use_mibig:
        name_to_idx = {}
        for cluster_idx, cluster_name in enumerate(bgc_collection.bgc_name_tuple):
            name_to_idx[cluster_name] = cluster_idx

        for bgc in mibig_set:
            mibig_set_indices.add(name_to_idx[bgc])
    return mibig_set_indices
