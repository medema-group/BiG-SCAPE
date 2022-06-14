import logging
import os

def generate_pfam_colors_matrix(pfam_domain_colors):
    '''
    :param pfam_domain_colors: tab-delimited file
    :return: dictionary with pfam ID as key and rgb colors as value
    '''
    pfam_colors = {}

    if os.path.isfile(pfam_domain_colors):
        print("  Found file with Pfam domain colors")
        with open(pfam_domain_colors, "r") as cat_handle:
            for line in cat_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    domain = row[0]
                    rgb = row[-1]
                    pfam_colors[domain] = rgb
    else:
        logging.warning("  File pfam_domain_colors was NOT found")

    return pfam_colors
