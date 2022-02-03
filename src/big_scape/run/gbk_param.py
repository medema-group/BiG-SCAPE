"""Module containing GBK parameter helper functions and classes

Author: Arjan Draisma
"""
import logging
import sys

class GbkParam:
    """Class which keeps track of run options relating to Gbk files"""

    include: list
    exclude: list

    def __init__(self, options):
        # set included gbk files
        self.set_include_list(options)

        # set excluded gbk files
        self.set_exclude_list(options)

    def set_include_list(self, options):
        """Sets the gbk inclusion list based on given options

        Inputs:
        - options: options object from CMD_parser"""

        # Exclude single string
        self.include = options.include_gbk_str
        if len(self.include) == 1 and self.include[0] == "*":
            logging.info(" Including all files")
        elif len(self.include) == 1 and self.include[0] == "":
            logging.error(" no strings specified for '--include_gbk_str'")
            sys.exit(1)
        else:
            logging.info(" Including files with one or more of the following strings in their filename: \
                '{}'".format("', '".join(self.include)))

    def set_exclude_list(self, options):
        """Sets the gbk exclusion list based on given options

        Inputs:
        - options: options object from CMD_parser"""

        self.exclude = options.exclude_gbk_str
        if self.exclude != []:
            logging.info(" Skipping files with one or more of the following strings in their filename: \
                '{}'".format("', '".join(self.exclude)))
