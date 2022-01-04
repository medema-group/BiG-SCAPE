"""Module containing network parameter helper functions and classes

Author: Arjan Draisma
"""

import os
import src.utility

class NetworkParam():
    """Class which keeps track of run options relating to generation of networks"""

    # anchor domains for networking
    anchor_domains: set

    def __init__(self, options):
        self.set_anchor_domains(options.anchorfile)

    def set_anchor_domains(self, anchorfile: str):
        """Set anchor domain data from file"""
        if os.path.isfile(anchorfile):
            self.anchor_domains = src.utility.get_anchor_domains(anchorfile)
        else:
            self.anchor_domains = set()
            # TODO: convert to log warning
            print("File with list of anchor domains not found")
