"""Module containing io parameter helper functions and classes

Author: Arjan Draisma
"""

import sys
import os

import src.utility as utility

class DirParam():
    """Class which keeps track of run options relating to directories and cache locations"""
    input: str
    output: str
    pfam: str
    mibig: str

    cache: str
    bgc_fasta: str
    domtable: str
    pfs: str
    pfd: str
    domains: str

    def __init__(self, options):
        # input dir
        self.set_input_dir(options)

        # output dir
        self.set_output_dir(options)
        self.prepare_output_dir()

        # cache dir
        self.set_cache_dir(options)
        self.prepare_cache_dir()

        # log dir
        self.set_log_dir(options)
        self.prepare_log_dir()

        # pfam dir
        self.set_pfam_dir(options)

        return

    def set_input_dir(self, options):
        """Checks the structure of input folder, and checks if valid data exists

        Inputs:
        - options: options object from CMD_parser"""

        # TODO: check if data exists
        self.input = options.inputdir

    def set_output_dir(self, options):
        """Checks if an output folder was given, and checks if it is writeable

        Inputs:
        - options: options object from CMD_parser"""

        if options.outputdir == "":
            # TODO: convert to log error
            print("please provide a name for an output folder using parameter -o or --outputdir")
            sys.exit(0)
        # TODO: check if writable
        self.output = options.outputdir

    def set_pfam_dir(self, options):
        """Checks if all necessary Pfam files exist in Pfam folder

        Inputs:
        - options: options object from CMD_parser"""

        h3f_exists = os.path.isfile(os.path.join(options.pfam_dir, "Pfam-A.hmm.h3f"))
        h3i_exists = os.path.isfile(os.path.join(options.pfam_dir, "Pfam-A.hmm.h3i"))
        h3m_exists = os.path.isfile(os.path.join(options.pfam_dir, "Pfam-A.hmm.h3m"))
        h3p_exists = os.path.isfile(os.path.join(options.pfam_dir, "Pfam-A.hmm.h3p"))
        if not (h3f_exists and h3i_exists and h3m_exists and h3p_exists):
            print("One or more of the necessary Pfam files (.h3f, .h3i, .h3m, .h3p) \
                were not found")
            if os.path.isfile(os.path.join(options.pfam_dir, "Pfam-A.hmm")):
                print("Please use hmmpress with Pfam-A.hmm")
            else:
                print("Please download the latest Pfam-A.hmm file from http://pfam.xfam.org/")
                print("Then use hmmpress on it, and use the --pfam_dir parameter \
                    to point to the location of the files")
            sys.exit()
        else:
            self.pfam = options.pfam_dir

    def set_cache_dir(self, options):
        """Sets cache folder associated with this run
        Creates new folders if necessary

        Inputs:
        - options: options object from CMD_parser"""

        self.cache = os.path.join(self.output, "cache")
        self.bgc_fasta = os.path.join(self.cache, "fasta")
        self.domtable = os.path.join(self.cache, "domtable")
        self.pfs = os.path.join(self.cache, "pfs")
        self.pfd = os.path.join(self.cache, "pfd")
        self.domains = os.path.join(self.cache, "domains")


    def set_log_dir(self, options):
        """Sets log directory associated with this run
        Creates a new directory if necessary

        Inputs:
        - options: options object from CMD_parser"""
        self.log = os.path.join(self.output, "logs")

    def prepare_output_dir(self):
        """Prepares the output directory by creating new folders"""
        # create output directory within output directory
        # TODO: simplify
        utility.create_directory(self.output, "Output", False)

    def prepare_cache_dir(self):
        """Prepares the output directory by creating new folders"""
        # create output directory within output directory
        # TODO: simplify
        utility.create_directory(self.cache, "Cache", False)
        utility.create_directory(self.bgc_fasta, "BGC fastas", False)
        utility.create_directory(self.domtable, "Domtable", False)
        utility.create_directory(self.domains, "Domains", False)
        utility.create_directory(self.pfs, "pfs", False)
        utility.create_directory(self.pfd, "pfd", False)

    def prepare_log_dir(self):
        """Prepares the output directory by creating new folders"""
        # TODO: remove?
        utility.create_directory(self.log, "Logs", False)
        # TODO: move elsewhere
        utility.write_parameters(self.log, sys.argv)

