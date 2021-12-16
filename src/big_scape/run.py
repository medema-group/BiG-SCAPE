"""Module containing classes to keep track of current BiG-SCAPE run details

Author: Arjan Draisma
"""
import sys
import os
from time import struct_time
import time

from src.utility.misc import get_anchor_domains

class Run:
    """Class to keep track of important run-specific details, based on given options
    """
    # I/O
    input_dir: str
    output_dir: str
    pfam_dir: str

    # anchor domain set
    anchor_domains: set = set()

    # mibig
    use_mibig = False

    # run mode
    run_mode: str

    # for logging
    start_time: struct_time
    run_name: str
    run_label: str
    run_data: dict

    def __init__(self, options):
        # input dir
        self.set_input_dir(options)

        # output dir
        self.set_output_dir(options)

        # pfam dir
        self.set_pfam_files(options)

        # mibig
        self.set_mibig(options)

        # run mode
        self.set_run_mode(options)

    def get_cache_folder(self):
        """Returns cache folder associated with this run

        returns: path to folder for cache files
        """
        return self.output_dir + "/cache"

    def set_anchor_domains(self, anchorfile: str):
        """Set anchor domain data from file"""
        if os.path.isfile(anchorfile):
            self.anchor_domains = get_anchor_domains(anchorfile)
        else:
            # TODO: convert to log warning
            print("File with list of anchor domains not found")

    def set_input_dir(self, options):
        """Checks the structure of input folder, and checks if valid data exists

        Inputs:
        - options: options object from CMD_parser"""
        # TODO: check if data exists
        self.input_dir = options.inputdir

    def set_output_dir(self, options):
        """Checks if an output folder was given, and checks if it is writeable

        Inputs:
        - options: options object from CMD_parser"""

        if options.outputdir == "":
            # TODO: convert to log error
            print("please provide a name for an output folder using parameter -o or --outputdir")
            sys.exit(0)
        # TODO: check if writable
        self.output_dir = options.outputdir

    def set_pfam_files(self, options):
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
            self.pfam_dir = options.pfam_dir

    def set_mibig(self, options):
        """Checks if the Mibig options that were passed are valid, and if so,
        returns if mibig is to be used

        Inputs:
        - options: options object from CMD_parser"""

        selected_mibig = 0
        if options.mibig21:
            selected_mibig += 1
        if options.mibig14:
            selected_mibig += 1
        if options.mibig13:
            selected_mibig += 1

        if selected_mibig > 1:
            sys.exit("Error: choose only one MIBiG version")

        self.use_mibig = selected_mibig == 1

    def set_run_mode(self, options):
        """Parses and sets the run mode of this run from options

        Inputs:
        - options: options object from CMD_parser"""
        run_mode_string = ""
        networks_folder_all = "networks_all"
        if options.hybrids:
            networks_folder_all += "_hybrids"
            run_mode_string += "_hybrids"
        if options.mode == "auto":
            networks_folder_all += "_auto"
            run_mode_string += "_auto"
        elif options.mode == "glocal":
            networks_folder_all += "_glocal"
            run_mode_string += "_glocal"
        else:
            run_mode_string += "_global"

        self.run_mode = run_mode_string


    def start(self, options):
        """Start the run: set a run name and record the start time"""
        # start time
        self.start_time = time.time()

        localtime = time.localtime(self.start_time)
        # generate run name
        self.run_name = "{}{}".format(time.strftime("%Y-%m-%d_%H-%M-%S", localtime), self.run_mode)

        if options.label:
            self.run_name = self.run_name + "_" + options.label

        # record run data
        # TODO: find out whether this is needed in this way
        self.run_data = {}
        self.run_data["start_time"] = time.strftime("%d/%m/%Y %H:%M:%S", localtime)
        self.run_data["parameters"] = " ".join(sys.argv[1:])
        self.run_data["input"] = {}
