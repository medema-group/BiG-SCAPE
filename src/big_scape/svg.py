import os

from glob import glob

import src.utility as utility




def get_available_svg(run):
    available_svg = set()
    for svg in glob(os.path.join(run.directories.svg, "*.svg")):
        (root, ext) = os.path.splitext(svg)
        available_svg.add(root.split(os.sep)[-1])
    return available_svg

def generate_images(run, cluster_base_names, gen_bank_dic, pfam_info, bgc_info):
    # All available SVG files
    available_svg = get_available_svg(run)

    # Which files actually need to be generated
    working_set = cluster_base_names - available_svg

    if len(working_set) > 0:
        color_genes = {}
        color_domains = utility.read_color_domains_file()
        pfam_domain_categories = {}

        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        print("  Reading BGC information and writing SVG")
        for bgc in working_set:
            with open(gen_bank_dic[bgc][0], "r") as handle:
                utility.SVG(False, os.path.join(run.directories.svg, bgc+".svg"), handle, bgc,
                            os.path.join(run.directories.pfd, bgc+".pfd"), True, color_genes,
                            color_domains, pfam_domain_categories, pfam_info,
                            bgc_info[bgc].records, bgc_info[bgc].max_width)

        color_genes.clear()
        color_domains.clear()
        pfam_domain_categories.clear()
    elif len(working_set) == 0:
        print("  All SVG from the input files seem to be in the SVG folder")

    available_svg.clear()
