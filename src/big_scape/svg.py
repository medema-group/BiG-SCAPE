import logging
import os

from glob import glob
from src.data.database import Database
from src.data import get_cluster_gbk_dict, gen_bgc_info_for_svg, get_cluster_name_list

import src.utility as utility


def get_available_svg(run):
    """Gets a set of available svg files from the svg directory"""
    available_svg = set()
    for svg in glob(os.path.join(run.directories.svg, "*.svg")):
        (root, ext) = os.path.splitext(svg)
        available_svg.add(root.split(os.sep)[-1])
    return available_svg

def generate_images(run, database: Database, pfam_info):
    """Generates SVG images
    
    Inputs:
        run: run details for this execution of BiG-SCAPE
        cluster_base_names: base names of all clusters in this run
        gen_bank_dict: dictionary of gen bank gbk files used
        pfam_info: pfam info from parsing Pfam-A.hmm
        bgc_info: dictionary of format {bgcname: {field: obj}}
    """
    gen_bank_dict = get_cluster_gbk_dict(run, database)

    bgc_info = gen_bgc_info_for_svg(database)

    bgc_names = get_cluster_name_list(database)

    # All available SVG files
    available_svg = get_available_svg(run)

    # Which files actually need to be generated
    working_set = set(bgc_names) - available_svg

    if len(working_set) > 0:
        color_genes = {}
        color_domains = utility.read_color_domains_file()
        pfam_domain_categories = {}

        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        logging.info("  Reading BGC information and writing SVG")
        for bgc in working_set:
            with open(gen_bank_dict[bgc], "r") as handle:
                utility.SVG(False, os.path.join(run.directories.svg, bgc+".svg"), handle, bgc,
                            database, color_genes,
                            color_domains, pfam_domain_categories, pfam_info,
                            bgc_info[bgc]["records"], bgc_info[bgc]["max_width"])

        color_genes.clear()
        color_domains.clear()
        pfam_domain_categories.clear()
    elif len(working_set) == 0:
        logging.info("  All SVG from the input files seem to be in the SVG folder")

    available_svg.clear()
