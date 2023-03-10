"""Module containing code to load and store GBK files"""

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from src.genbank.region import Region


class GBK:
    """Class to describe a given GBK file"""

    path: Path
    metadata: dict[str, str]
    region: Region

    def __init__(self, path):
        self.path = path
        self.metadata = {}

    @classmethod
    def parse(cls, path: Path):
        """Parses a GBK file and returns a GBK object with all necessary information"""
        gbk = cls(path)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))

        # go through features
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                region = Region.create_region(feature)
                gbk.region = region

        return gbk
