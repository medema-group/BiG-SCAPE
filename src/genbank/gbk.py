"""Module containing code to load and store GBK files"""

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError
from src.genbank.region import Region


class GBK:
    """
    Class to describe a given GBK file

    Attributes:
        path: Path
        metadata: Dict[str, str]
        region: Region
    """

    def __init__(self, path):
        self.path = path
        self.metadata = {}
        self.region: Region = None

    @classmethod
    def parse(cls, path: Path):
        """Parses a GBK file and returns a GBK object with all necessary information"""
        gbk = cls(path)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))

        # go through features
        # TODO load, index, populate objects (load into dicts, loop dict and populate)
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                if gbk.region is not None:
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(feature)
                gbk.region = region

        return gbk
