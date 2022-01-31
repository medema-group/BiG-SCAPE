import math
from src.big_scape.bgc_info import BgcInfo

def generate_domain_name_info(bgc: BgcInfo):
    bgc.domain_name_info = {}
    for cluster_name in bgc.ordered_domain_list:
        for i in range(10):
            strand = "+" if i % 2 == 0 else "-"
            bgc.domain_name_info[cluster_name] = f"{bgc.name}_ORF1:gid::pid::loc:{i*10}:{i+1*10}:strand:{strand}"
    return bgc

def add_similar_bgc_domains(bgc_a, bgc_b, num_domains):
    # 10% of domains won't correspond (tail ends)
    margins = math.ceil(num_domains / 10)
    # we assign from this total number of domains
    total_domains = num_domains * 2 + margins
    for i in range(total_domains):
        dom_name = "DOM_" + str(i)
        # bottom half are in a
        if i < total_domains - margins / 2:
            bgc_a.ordered_domain_list.append(dom_name)
        # top half are in b
        if i > margins / 2:
            bgc_b.ordered_domain_list.append(dom_name)
    return bgc_a, bgc_b


def add_distinct_bgc_domains(bgc_a, bgc_b, num_domains):
    """This function adds domains to each bgc, and ensures none of the domains are shared
    
    This creates two highly dissimilar gene clusters
    
    Input:
        bgc_a, bgc_b: BgcInfo - gene clusters to fill with domains
        num_domains: int - Total number of domains to add to each bgc
    """
    # total amount of domains is requested domains * 2
    # even numbers are assigned to A, uneven numbers are assigned to B
    for i in range(num_domains*2):
        dom_name = "DOM_" + str(i)
        # even are in a
        if i % 2 == 0:
            bgc_a.ordered_domain_list.append(dom_name)
        # odd are in b
        else:
            bgc_b.ordered_domain_list.append(dom_name)
    return bgc_a, bgc_b

def create_cluster_couple(similar: bool, num_domains: int):
    bgc_a = BgcInfo("test_similar_bgc_a")
    bgc_b = BgcInfo("test_similar_bgc_b")

    bgc_a.ordered_domain_list = []
    bgc_b.ordered_domain_list = []

    if similar:
        add_similar_bgc_domains(bgc_a, bgc_b, num_domains)
    else:
        add_distinct_bgc_domains(bgc_a, bgc_b, num_domains)

    # init sets
    bgc_a.ordered_domain_set = set(bgc_a.ordered_domain_list)
    bgc_b.ordered_domain_set = set(bgc_b.ordered_domain_list)

    generate_domain_name_info(bgc_a)
    generate_domain_name_info(bgc_b)

    return bgc_a, bgc_b