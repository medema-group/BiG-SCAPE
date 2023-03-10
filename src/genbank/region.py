"""Module containing code to load and store AntiSMASH regions"""

import logging

from Bio.SeqFeature import SeqFeature


class Region:
    """Class to describe a region within an Antismash GBK"""

    number: int

    def __init__(self, number):
        self.number = number

    @classmethod
    def create_region(cls, feature: SeqFeature):
        """Creates a region object from a region feature in a GBK file"""
        if feature.type != "region":
            logging.error(
                "Feature is not of correct type! (expected: region, was: %s)",
                feature.type,
            )
            raise ValueError()

        if "region_number" not in feature.qualifiers:
            logging.error("region_number qualifier not found in region feature!")
            raise ValueError()

        region_number = int(feature.qualifiers["region_number"][0])

        region = cls(region_number)

        return region
