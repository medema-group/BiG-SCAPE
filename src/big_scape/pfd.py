import os
from array import array
from collections import defaultdict

def parse_pfd(run, cluster_base_names, gene_domain_count, corebiosynthetic_pos,
              bgc_gene_orientation, bgc_info):
    pfd_dict_domains = defaultdict(int)
    orf_keys = {}
    for outputbase in cluster_base_names:
        gene_domain_count[outputbase] = array('B')
        corebiosynthetic_pos[outputbase] = array('H')
        bgc_gene_orientation[outputbase] = array('b')
        pfd_file = os.path.join(run.directories.pfd, outputbase + ".pfd")

        #pfd_dict_domains contains the number of domains annotated in the
        # pfd file for each orf tag
        with open(pfd_file, "r") as pfdf:
            for line in pfdf:
                pfd_dict_domains[line.strip().split("\t")[-1]] += 1

        # extract the orf number from the tag and use it to traverse the BGC
        for orf in pfd_dict_domains.keys():
            orf_num = int(orf.split(":")[0].split("_ORF")[1])
            orf_keys[orf_num] = orf

        orf_num = 0
        for orf_key in sorted(orf_keys.keys()):
            orf = orf_keys[orf_key]
            if orf[-1] == "+":
                bgc_gene_orientation[outputbase].append(1)
            else:
                bgc_gene_orientation[outputbase].append(-1)

            gene_domain_count[outputbase].append(pfd_dict_domains[orf])

            if orf in bgc_info[outputbase].biosynthetic_genes:
                corebiosynthetic_pos[outputbase].append(orf_num)
            orf_num += 1

        pfd_dict_domains.clear()
        orf_keys.clear()

        ## TODO: if len(corebiosynthetic_position[outputbase]) == 0
        ## do something with the list of pfam ids. Specifically, mark
        ## (in this case TODO or always?) as biosynthetic genes, the ones that contain
        ## domains from a special list. This list of special domains
        ## comes from predicted domains within the CDSs marked as 'sec_met'
        ## by antismash
