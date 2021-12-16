"""Module containing io parameter helper functions and classes

Author: Arjan Draisma
"""

import sys
import os

class DirParam():
    input: str
    output: str
    cache: str
    pfam: str
    mibig: str

    def __init__(self, options):
        # input dir
        self.set_input_dir(options)

        # output dir
        self.set_output_dir(options)

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

    def set_cache_folder(self, options):
        """Sets cache folder associated with this run
        Creates a new folder if necessary
        
        Inputs:
        - options: options object from CMD_parser"""
        self.cache = self.output + "/cache"
