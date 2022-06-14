"""Module containing mibig parameter helper functions and classes

Author: Arjan Draisma
"""

import logging
import sys
import os

class MibigParam():
    """Class which keeps track of run options relating to mibig usage"""
    # whether to use mibig in a run
    use_mibig = False

    file_name: str
    expected_num_bgc: int

    zip_path: str
    gbk_path: str

    def __init__(self, options):
        self.set_mibig(options)
        self.set_mibig_file_path(options)

    def set_mibig(self, options):
        """Checks if the Mibig options that were passed are valid, and if so,
        sets the use mibig flag if it is to be used

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
            logging.error("choose only one MIBiG version")
            sys.exit(1)

        self.use_mibig = selected_mibig == 1


    def set_mibig_file_path(self, options):
        """Sets the file name base and expected MiBIG BGC count if MiBIG is
        used

        Inputs:
        - options: options object from CMD_parser"""

        if not self.use_mibig:
            return

        if options.mibig21:
            self.file_name = "MIBiG_2.1_final"
            self.expected_num_bgc = 1923
        elif options.mibig14:
            self.file_name = "MIBiG_1.4_final"
            self.expected_num_bgc = 1808
        else:
            self.file_name = "MIBiG_1.3_final"
            self.expected_num_bgc = 1393

        self.zip_path = os.path.join(
            options.mibig_path, self.file_name + ".zip"
        )
        self.gbk_path = os.path.join(options.mibig_path, self.file_name)
