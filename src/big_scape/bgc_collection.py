from typing import Dict

import logging

from src.big_scape.bgc_info import BgcInfo



class BgcCollection:
    """Class to contain a collection of BGC information"""
    # a list of all bgcs recorded in this collection
    bgc_name_list: list
    # the same list as a set
    bgc_name_set: set
    # a tuple needed later
    # TODO: refactor out
    bgc_name_tuple: tuple

    # dictionary of domains per bgc
    bgc_ordered_domain_list: dict

    # the dictionary of actual bgc info objects
    bgc_collection_dict: Dict[str, BgcInfo]

    def initialize(self, cluster_name_list):
        """Initialize the collection

        Input:
            cluster_name_list: a list of names of all clusters to be added to this collection
        """
        self.bgc_name_list = cluster_name_list
        self.bgc_name_set = set(cluster_name_list)
        self.bgc_name_tuple = tuple(sorted(cluster_name_list))

        self.bgc_collection_dict = {}
        for cluster_name in cluster_name_list:
            self.bgc_collection_dict[cluster_name] = BgcInfo(cluster_name)

    def add_bgc_info(self, bgc_info_dict):
        """Add BGC information to each cluster in this collection

        Inputs:
            bgc_info_dict: bgc info dictionary from gbk.import_gbks
        """
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in bgc_info_dict:
                bgc_info = bgc_info_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].bgc_info = bgc_info
            else:
                # TODO: replace with logger
                logging.warning("BGC info for %s was not found", cluster_name)

    def add_source_gbk_files(self, source_gbk_file_dict):
        """Add source GBK files to the collection

        Inputs:
            source_gbk_file_dict: source file dict from gbk.import_gbks"""
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in source_gbk_file_dict:
                source_gbk_file = source_gbk_file_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].src_gbk_file = source_gbk_file
            else:
                # TODO: replace with logger
                logging.warning("BGC info for %s was not found", cluster_name)

    def add_gene_domain_counts(self, gene_domain_count_dict):
        """Add gene domain counts to bgcs

        Inputs:
            gene_domain_count_dict: {bgc: [domain counts]} dictionary from big_scape.parse_pfd
        """
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in gene_domain_count_dict:
                gene_domain_counts = gene_domain_count_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].num_genes = len(gene_domain_counts)
                self.bgc_collection_dict[cluster_name].gene_domain_counts = gene_domain_counts
            else:
                self.bgc_collection_dict[cluster_name].num_genes = 0
                # TODO: replace with logger
                logging.warning("Domain count info for %s was not found", cluster_name)

    def add_gene_orientations(self, gene_orientation_dict):
        """Add gene orientation information to bgcs

        Inputs:
            gene_orientation_dict: {bgc: [orientations]} dictionary from big_scape.parse_pfd"""
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in gene_orientation_dict:
                gene_orientations = gene_orientation_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].gene_orientations = gene_orientations
            else:
                # TODO: replace with logger
                logging.warning("Gene orientation info for %s was not found", cluster_name)

    def add_bio_synth_core_pos(self, bio_synth_core_positions):
        """Add positions of core biosynthetic genes for each bgc to this collection
        
        Inputs:
            bio_synth_core_positions: {bgc: [positions]} dictionary from big_scape.parse_pfd
        """
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in bio_synth_core_positions:
                cluster_bio_synth_core_positions = bio_synth_core_positions[cluster_name]
                self.bgc_collection_dict[cluster_name].bio_synth_core_positions = cluster_bio_synth_core_positions
            else:
                # TODO: replace with logger
                logging.warning("Biosynthetic gene core position info for %s was not found", cluster_name)

    def init_gene_strings(self):
        for cluster_name in self.bgc_collection_dict:
            self.bgc_collection_dict[cluster_name].init_gene_string()
