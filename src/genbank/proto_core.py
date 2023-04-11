"""Module containing code to load and store AntiSMASH protoclusters"""

# from python
import logging
from typing import Any, Optional

# from dependencies
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError

# from this module
from src.genbank.bgc_record import BGCRecord


class ProtoCore(BGCRecord):
    """
    Class to describe a protocore within an Antismash GBK

    Attributes:
        nt_start: int
        nt_stop: int
        product: str
        number: int
    """

    def __init__(self, number: int):
        super().__init__()
        self.number = number

    def save(self, commit=True):
        """Stores this protocore in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        return super().save("proto_core", commit)

    @classmethod
    def parse(cls, feature: SeqFeature, parent_gbk: Optional[Any] = None):
        """Creates a Protocore object from a region feature in a GBK file

        Args:
            feature (SeqFeature): proto_core antiSMASH genbank feature

        Raises:
            InvalidGBKError: invalid or missing field

        Returns:
            ProtoCore: protocore object
        """
        if feature.type != "proto_core":
            logging.error(
                "Feature is not of correct type! (expected: proto_core, was: %s)",
                feature.type,
            )
            raise InvalidGBKError()

        if "protocluster_number" not in feature.qualifiers:
            logging.error(
                "protocluster_number qualifier not found in proto_core feature!"
            )
            raise InvalidGBKError()

        proto_core_number = int(feature.qualifiers["protocluster_number"][0])

        proto_core = cls(proto_core_number)
        proto_core.parse_bgc_record(feature, parent_gbk=parent_gbk)

        return proto_core
