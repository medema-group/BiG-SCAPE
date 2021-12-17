"""Module containing cluster parameter helper functions and classes

Author: Arjan Draisma
"""

import sys

class ClusterParam():
    """Class which keeps track of run options relating to clustering"""
    cutoff_list: list
    max_cutoff: int

    def __init__(self, options):
        # add regular cutoffs
        self.add_cutoffs(options)

        # if we want to classify by clans make sure that the clanCutoff is included
        # in the cutoffs to do AP in clusterJsonBatch
        if options.clans:
            self.add_clan_cutoffs(options)

    def add_cutoffs(self, options):
        self.cutoff_list = options.cutoffs
        for cutoff in self.cutoff_list:
            if cutoff <= 0.0 or cutoff > 1.0:
                sys.exit(" Invalid cutoff value {}".format(str(cutoff)))
        self.max_cutoff = max(self.cutoff_list)

    def add_clan_cutoffs(self, options):
        family_cutoff, clan_cutoff = options.clan_cutoff
        if clan_cutoff < family_cutoff:
            sys.exit("Error: first value in the clan_cutoff parameter should be smaller than the \
            second")
        if family_cutoff not in self.cutoff_list:
            if family_cutoff <= 0.0 or family_cutoff > 1.0:
                sys.exit("Error: invalid cutoff value for GCF calling")
            else:
                self.cutoff_list.append(family_cutoff)
                self.cutoff_list.sort()

        if clan_cutoff <= 0.0 or clan_cutoff > 1.0:
            sys.exit("Error: invalid cutoff value for GCC calling")
