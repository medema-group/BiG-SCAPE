"""Module containing classes to keep track of current BiG-SCAPE run details

Author: Arjan Draisma
"""
import logging
import sys
import os
import time

from .dir_param import DirParam
from .gbk_param import GbkParam
from .mibig_param import MibigParam
from .distance_param import DistParam
from .cluster_param import ClusterParam
from .network_param import NetworkParam

class Run:
    """Class to keep track of important run-specific details, based on given
    options
    """
    ## subsections

    # options used to generate this run
    options: object

    # files
    directories: DirParam
    gbk: GbkParam

    # mibig
    mibig: MibigParam

    # distance calculation
    distance: DistParam

    # clustering
    cluster: ClusterParam

    # networking
    network: NetworkParam

    ## other run parameters
    # run mode
    run_mode: str

    # domain include list
    has_includelist: bool
    domain_includelist: set

    # valid/banned classes
    valid_classes: set
    user_banned_classes: set

    # for logging
    start_time: time.struct_time
    run_name: str
    run_label: str
    run_data: dict

    def set_run_mode(self, options):
        """Parses and sets the run mode of this run from options

        Inputs:
        - Options: options object from cmd_parser
        """
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

    def set_domain_includelist(self, options):
        """Sets flag to use includelist to true if it is present in options

        Inputs:
        - Options: options object from cmd_parser
        """
        self.has_includelist = False
        if options.domain_includelist:
            bigscape_path = os.path.dirname(os.path.realpath(__file__))
            includelist_path = os.path.join(
                bigscape_path,
                "domain_includelist.txt"
            )
            if os.path.isfile(includelist_path):
                self.domain_includelist = set()
                # file to read includelist from
                for line in open(includelist_path, "r", encoding="utf-8"):
                    if line[0] == "#":
                        continue
                    self.domain_includelist.add(line.split("\t")[0].strip())
                if len(self.domain_includelist) == 0:
                    log_line = ("--domain_includelist used, but no "
                    "domains found in the file")
                    logging.warning(log_line)
                else:
                    self.has_includelist = True
            else:
                logging.error("domain_includelist.txt file not found")
                sys.exit(1)

    def set_valid_classes(self, options):
        """Function to set which classes will be analyzed (if in the
        options_classify mode)

        Inputs:
        - Options: options object from cmd_parser
        """
        self.valid_classes = set()
        for key in self.distance.bgc_class_weight:
            self.valid_classes.add(key.lower())
        self.user_banned_classes = set(
            [a.strip().lower() for a in options.banned_classes]
        )
        self.valid_classes = self.valid_classes - self.user_banned_classes


    def init(self, options):
        """Function to initialize run parameters based on given options object

        Inputs:
        - Options: options object from cmd_parser
        """
        self.options = options

        self.directories = DirParam(options)
        self.gbk = GbkParam(options)

        self.mibig = MibigParam(options)
        self.distance = DistParam(options)
        self.cluster = ClusterParam(options)
        self.network = NetworkParam(options)

        self.set_run_mode(options)

        self.set_domain_includelist(options)

        self.set_valid_classes(options)

    def start(self, skip_dir=False):
        """Start the run: set a run name and record the start time

        Inputs:
        - options: options object from CMD_parser"""
        # start time
        self.start_time = time.time()

        localtime = time.localtime(self.start_time)
        # generate run name
        timestamp = time.strftime("%Y-%m-%d_%H-%M-%S", localtime)
        self.run_name = f"{timestamp}{self.run_mode}"

        if self.options.label:
            self.run_name = self.run_name + "_" + self.options.label

        # record run data
        self.run_data = {}
        self.run_data["start_time"] = time.strftime(
            "%d/%m/%Y %H:%M:%S",
            localtime
        )
        self.run_data["parameters"] = " ".join(sys.argv[1:])
        self.run_data["input"] = {}

        if skip_dir:
            return

        self.directories.set_run_dependent_dir(self.run_name)
        self.directories.prepare_run_dependent_dir()

    def end(self):
        """Call and store the end time of the run"""
        end_time = time.time()
        duration = end_time - self.start_time
        self.run_data["end_time"] = time.strftime(
            "%d/%m/%Y %H:%M:%S",
            time.localtime(end_time)
        )
        hours = (duration // 3600)
        minutes = ((duration % 3600) // 60)
        seconds = ((duration % 3600) % 60)
        self.run_data["duration"] = f"{hours}h{minutes}m{seconds}s"

    def report_runtime(self):
        """Report the runtime of this run to the logger"""

        runtime = time.time()-self.start_time
        runtime_string = f"Main function took {runtime:.3f} s"

        runtime_path = os.path.join(self.directories.log, "runtimes.txt")

        with open(runtime_path, 'a', encoding="utf-8") as timings_file:
            timings_file.write(runtime_string + "\n")

        # print runtime
        logging.info(runtime_string)
