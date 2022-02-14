import os
import sys

from typing import Dict

import debugpy
import logging
import pickle

from src.big_scape.util import get_ordered_domain_list
from src.big_scape.bgc_info import BgcInfo
from src.bgctools import generate_domain_name_info_dict
from src.utility import fasta_parser, save_domain_seqs



class BgcCollection:
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
        self.bgc_name_list = cluster_name_list
        self.bgc_name_set = set(cluster_name_list)
        self.bgc_name_tuple = tuple(sorted(cluster_name_list))

        self.bgc_collection_dict = {}
        for cluster_name in cluster_name_list:
            self.bgc_collection_dict[cluster_name] = BgcInfo(cluster_name)

    def init_ordered_domain_list(self, run):
        # Get the ordered list of domains
        self.bgc_ordered_domain_list = get_ordered_domain_list(run, self.bgc_name_set)
        for cluster_name in self.bgc_ordered_domain_list:
            domain_list = self.bgc_ordered_domain_list[cluster_name]
            self.bgc_collection_dict[cluster_name].ordered_domain_list = domain_list
        
            # we will also want a set later on for various uses
            self.bgc_collection_dict[cluster_name].ordered_domain_set = set(domain_list)

    def add_bgc_info(self, bgc_info_dict):
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in bgc_info_dict:
                bgc_info = bgc_info_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].bgc_info = bgc_info
            else:
                # TODO: replace with logger
                logging.warning("BGC info for %s was not found", cluster_name)

    def add_source_gbk_files(self, source_gbk_file_dict):
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in source_gbk_file_dict:
                source_gbk_file = source_gbk_file_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].src_gbk_file = source_gbk_file
            else:
                # TODO: replace with logger
                logging.warning("BGC info for %s was not found", cluster_name)

    def add_gene_domain_counts(self, gene_domain_count_dict):
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
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in gene_orientation_dict:
                gene_orientations = gene_orientation_dict[cluster_name]
                self.bgc_collection_dict[cluster_name].gene_orientations = gene_orientations
            else:
                # TODO: replace with logger
                logging.warning("Gene orientation info for %s was not found", cluster_name)

    def add_bio_synth_core_pos(self, bio_synth_core_positions):
        for cluster_name in self.bgc_collection_dict:
            if cluster_name in bio_synth_core_positions:
                cluster_bio_synth_core_positions = bio_synth_core_positions[cluster_name]
                self.bgc_collection_dict[cluster_name].bio_synth_core_positions = cluster_bio_synth_core_positions
            else:
                # TODO: replace with logger
                logging.warning("Biosynthetic gene core position info for %s was not found", cluster_name)

    def load_domain_names_from_dict(self, run):
        try:
            with open(os.path.join(run.directories.cache, "BGCs.dict"), "r") as bgc_file:
                domain_name_info_dict = pickle.load(bgc_file)
                for cluster_name in domain_name_info_dict:
                    domain_name_info = domain_name_info_dict[cluster_name]
                    self.bgc_collection_dict[cluster_name].domain_name_info = domain_name_info
                bgc_file.close()
        except IOError:
            logging.error("BGCs file not found...")
            sys.exit(1)

    def load_domain_names_from_pfd(self, run, try_resume_multiple_alignment):
        for outputbase in self.bgc_name_list:
            logging.debug("   Processing: %s", outputbase)

            pfd_file = os.path.join(run.directories.pfd, outputbase + ".pfd")
            filtered_matrix = [[part.strip() for part in l.split('\t')] for l in open(pfd_file)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(run.directories.bgc_fasta, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not try_resume_multiple_alignment:
                with open(fasta_file, "r") as fasta_file_handle:
                    # all fasta info from a BGC
                    fasta_dict = fasta_parser(fasta_file_handle)
                save_domain_seqs(filtered_matrix, fasta_dict,
                                         run.directories.domains, outputbase)

            domain_name_info = generate_domain_name_info_dict(filtered_matrix)
            self.bgc_collection_dict[outputbase].domain_name_info = domain_name_info

            del filtered_matrix[:]

    def save_domain_names_to_dict(self, run):
        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(run.directories.cache, "BGCs.dict"), "wb") as bgc_file:
            domain_name_info = {}
            for cluster_name in self.bgc_collection_dict:
                domain_name_info[cluster_name] = self.bgc_collection_dict[cluster_name].domain_name_info
            pickle.dump(domain_name_info, bgc_file)
            bgc_file.close()

    def init_gene_strings(self):
        for cluster_name in self.bgc_collection_dict:
            self.bgc_collection_dict[cluster_name].init_gene_string()

    def is_ready(self):
        return
