"""Module containing classes to keep track of current BiG-SCAPE run details

Author: Arjan Draisma
"""
import sys
import os
import time

from src.big_scape.run.dir_param import DirParam
from src.big_scape.run.mibig_param import MibigParam
from src.big_scape.run.pfam_param import PfamParam
from src.big_scape.run.cluster_param import ClusterParam
from src.big_scape.run.network_param import NetworkParam

class Run:
    """Class to keep track of important run-specific details, based on given options
    """
    ## subsections
    # io
    directories: DirParam

    # mibig
    mibig: MibigParam

    # pfam
    pfam: PfamParam

    # clustering
    cluster: ClusterParam

    # networking
    network: NetworkParam

    ## other run parameters
    # run mode
    run_mode: str
    has_query_bgc: bool

    # for logging
    start_time: time.struct_time
    run_name: str
    run_label: str
    run_data: dict

    def __init__(self, options):
        self.directories = DirParam(options)

        self.mibig = MibigParam(options)
        self.pfam = PfamParam(options)
        self.cluster = ClusterParam(options)
        self.network = NetworkParam(options)

        self.set_run_mode(options)

        self.set_has_query_bgc(options)


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

    def set_has_query_bgc(self, options):    
        self.has_query_bgc = False
        if options.query_bgc:
            self.has_query_bgc = True
            if not os.path.isfile(options.query_bgc):
                sys.exit("Error: Query BGC not found")


    def start(self, options):
        """Start the run: set a run name and record the start time

        Inputs:
        - options: options object from CMD_parser"""
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

