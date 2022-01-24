"""Module containing io parameter helper functions and classes

Author: Arjan Draisma
"""

import sys
import os

import src.utility as utility

class DirParam():
    """Class which keeps track of run options relating to directories and cache locations"""
    # folders
    input: str
    output: str
    pfam: str
    mibig: str

    # query bgc
    has_query_bgc: bool
    query_bgc: str
    query_bgc_name: str

    # cache
    cache: str
    bgc_fasta: str
    domtable: str
    pfs: str
    pfd: str
    domains: str

    # files and directories created per run
    # for network files
    network: str
    # network html
    network_html: str

    def __init__(self, options):
        # input dir
        self.set_input_dir(options)

        # bgc query folder
        self.set_has_query_bgc(options)

        # output dir
        self.set_output_dir(options)
        self.prepare_output_dir()

        # svg dir
        self.set_svg_dir()
        self.prepare_svg_dir()

        # cache dir
        self.set_cache_dir(options)
        self.prepare_cache_dir()

        # log dir
        self.set_log_dir(options)
        self.prepare_log_dir()

        # pfam dir
        self.set_pfam_dir(options)

    def set_input_dir(self, options):
        """Checks the structure of input folder, and checks if valid data exists

        Inputs:
        - options: options object from CMD_parser"""

        # TODO: check if data exists
        self.input = options.inputdir

    def set_has_query_bgc(self, options):
        """Set flag for query biosynthetic gene clusters & set query bgc folder

        Inputs:
        - options: options object from CMD_parser"""
        if options.query_bgc:
            self.has_query_bgc = True
            self.query_bgc_name = ".".join(options.query_bgc.split(os.sep)[-1].split(".")[:-1])
            if not os.path.isfile(options.query_bgc):
                sys.exit("Error: Query BGC not found")
            self.query_bgc = options.query_bgc
        else:
            self.has_query_bgc = False

    def set_output_dir(self, options):
        """Checks if an output folder was given, and checks if it is writeable
        Then sets output directory information on this class

        Inputs:
        - options: options object from CMD_parser"""

        if options.outputdir == "":
            # TODO: convert to log error
            print("please provide a name for an output folder using parameter -o or --outputdir")
            sys.exit(0)
        # TODO: check if writable
        self.output = options.outputdir

    def set_svg_dir(self):
        """Sets SVG output directory information on this class

        Inputs:
        - options: options object from CMD_parser"""

        # TODO: check if writable
        self.svg = os.path.join(self.output, "SVG")

    def set_run_dependent_dir(self, run_name):
        """Sets network output directory information on this class

        This is executed after Run.start()

        Inputs:
        - options: options object from CMD_parser"""

        # TODO: check if writable
        network = os.path.join(self.output, "network_files")
        # create output directory for network files *for this run*
        self.network = os.path.join(network, run_name)
        self.network_html = os.path.join(self.output, "html_content", "networks", run_name)

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
        """Prepares the output directory by creating new folder"""
        # TODO: simplify
        utility.create_directory(self.output, "Output", False)

    def prepare_svg_dir(self):
        """Prepares the svg directory by creating new folder"""
        utility.create_directory(self.svg, "SVG", False)

    def prepare_run_dependent_dir(self):
        """Prepares the network directory by creating new folder"""
        utility.create_directory(os.path.join(self.output, "network_files"), "Networks", False)
        utility.create_directory(self.network, "Network Files", False)

    def prepare_cache_dir(self):
        """Prepares the cache directory by creating new folders"""
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
