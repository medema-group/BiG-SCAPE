"""Module containing distance parameter helper functions and classes

Author: Arjan Draisma
"""

from collections import defaultdict

class DistParam:
    """Class which keeps track of run options relating to distance calculation"""
    bgc_class_weight: dict
    bgc_classes: dict
    bgc_class_names: tuple

    def __init__(self, options):
        self.set_weights(options)

    def set_weights(self, options):
        """Initializes weights with default options

        Parameters:
        - options: options object from CMD_parser
        """
        # TODO: implement custom options
        # Weights in the format J, DSS, AI, anchorboost
        # Generated with optimization results 2016-12-05.
        # Used the basic list of 4 anchor domains.
        self.bgc_class_weight = {}
        self.bgc_class_weight["PKSI"] = (0.22, 0.76, 0.02, 1.0)
        self.bgc_class_weight["PKSother"] = (0.0, 0.32, 0.68, 4.0)
        self.bgc_class_weight["NRPS"] = (0.0, 1.0, 0.0, 4.0)
        self.bgc_class_weight["RiPPs"] = (0.28, 0.71, 0.01, 1.0)
        self.bgc_class_weight["Saccharides"] = (0.0, 0.0, 1.0, 1.0)
        self.bgc_class_weight["Terpene"] = (0.2, 0.75, 0.05, 2.0)
        self.bgc_class_weight["PKS-NRP_Hybrids"] = (0.0, 0.78, 0.22, 1.0)
        self.bgc_class_weight["Others"] = (0.01, 0.97, 0.02, 4.0)

        # finally, define weights for mix

        # default when not separating in classes
        self.bgc_class_weight["mix"] = (0.2, 0.75, 0.05, 2.0)

        self.bgc_classes = defaultdict(list)

        # mix class will always be the last element of the tuple
        self.bgc_class_names = tuple(sorted(list(self.bgc_class_weight)) + ["mix"])
        assert self.bgc_class_names[-1] == 'mix'
