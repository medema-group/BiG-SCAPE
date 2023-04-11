"""Module containing code to load and store GBK files"""

# from python
import logging
from enum import Enum
from pathlib import Path
from typing import Dict, Optional, List


# from dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

# from other modules
from src.errors import InvalidGBKError
from src.data import DB

# from this module
from src.genbank.region import Region
from src.genbank.candidate_cluster import CandidateCluster
from src.genbank.proto_cluster import ProtoCluster
from src.genbank.proto_core import ProtoCore
from src.genbank.cds import CDS


class SOURCE_TYPE(Enum):
    QUERY = "query"
    MIBIG = "mibig"
    REFERENCE = "reference"


class GBK:
    """
    Class to describe a given GBK file

    Attributes:
        path: Path
        metadata: Dict[str, str]
        region: Region
        nt_seq: SeqRecord.seq
        genes: list[CDS]
        as_version: str
        source_type: SOURCE_TYPE
    """

    def __init__(self, path, source_type) -> None:
        self.path: Path = path
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None
        self.nt_seq: SeqRecord.seq = None
        self.genes: List[Optional[CDS]] = []
        self.as_version: Optional[str] = None
        self.source_type: SOURCE_TYPE = source_type

        # db specific fields
        self._db_id: Optional[int] = None

    def save(self, commit=True):
        """Stores this GBK in the database

        this returns the id of the GBK row

        Arguments:
            commit: commit immediately after executing the insert query"""
        gbk_table = DB.metadata.tables["gbk"]
        insert_query = (
            gbk_table.insert()
            .values(path=str(self.path), nt_seq=str(self.nt_seq))
            .returning(gbk_table.c.id)
            .compile()
        )

        # in the above query we add a returning statement. This makes it so that the
        # sqlite engine will be in the middle of a transaction (trying to give us the
        # returned value) and means we cannot commit. We will do this after digesting
        # the reutrn statement further on
        cursor_result = DB.execute(insert_query, False)

        # get return value
        return_row = cursor_result.fetchone()
        self._db_id = return_row[0]

        # only now that we have handled the return we can commit
        if commit:
            DB.commit()

    def save_all(self):
        """Stores this GBK and its children in the database. Does not commit immediately

        this function never commits"""
        self.save(False)
        self.region.save_all()

    @staticmethod
    def get_as_version(gbk_seq_record: SeqRecord):
        """Get AS version from GBK record

        Args:
            gbk_seq_record (SeqRecord): gbk seqrecord

        Returns:
            str: antismash version
        """

        try:
            as_version = gbk_seq_record.annotations["structured_comment"][
                "antiSMASH-Data"
            ]["Version"]
        except KeyError:
            # assume AS version 4
            as_version = "4"

        return as_version

    @classmethod
    def parse(cls, path: Path, source_type: SOURCE_TYPE):
        """Parses a GBK file and returns a GBK object with all necessary information

        Args:
            path (Path): path to gbk file

        Returns:
            GBK: GBK object
        """

        gbk = cls(path, source_type)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))
        gbk.nt_seq = record.seq

        as_version = GBK.get_as_version(record)
        gbk.as_version = as_version

        if int(as_version[0]) >= 5:
            gbk.parse_as5up(record)
        if int(as_version[0]) == 4:
            gbk.parse_as4(record)

        return gbk

    def parse_as4(self, record: SeqRecord):
        """Parses a GBK record of AS version 4 and returns a GBK object with all necessary information

        Args:
            record (SeqRecord): gbk file record

        Raises:
            InvalidGBKError: Invalid or missing fields in gbk record
        """

        # go through features, load into tmp dicts indexed by feature number
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "cluster":
                if self.region is not None:
                    # this should not happen, but just in case
                    # since we only have one region on an object
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(self, feature)
                self.region = region

            if feature.type == "CDS":
                cds = CDS.parse(feature)
                self.genes.append(cds)

    def parse_as5up(self, record: SeqRecord):
        """Parses a GBK record of AS versions 5 and up and returns a GBK object with all necessary information

        Args:
            record (SeqRecord): gbk file record

        Raises:
            InvalidGBKError: Invalid or missing fields in gbk record
        """

        tmp_cand_clusters = {}
        tmp_proto_clusters = {}
        tmp_proto_cores = {}

        # go through features, load into tmp dicts indexed by feature number
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                if self.region is not None:
                    # this should not happen, but just in case
                    # since we only have one region on an object
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(self, feature)
                self.region = region

            if feature.type == "cand_cluster":
                cand_cluster = CandidateCluster.parse(self, feature)
                tmp_cand_clusters[cand_cluster.number] = cand_cluster

            if feature.type == "protocluster":
                proto_cluster = ProtoCluster.parse(self, feature)
                tmp_proto_clusters[proto_cluster.number] = proto_cluster

            if feature.type == "proto_core":
                proto_core = ProtoCore.parse(self, feature)
                tmp_proto_cores[proto_core.number] = proto_core

            if feature.type == "CDS":
                cds = CDS.parse(feature)
                self.genes.append(cds)

        # add features to parent objects
        for proto_cluster_num, proto_cluster in tmp_proto_clusters.items():
            proto_cluster.add_proto_core(tmp_proto_cores[proto_cluster_num])

        for cand_cluster in tmp_cand_clusters.values():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                cand_cluster.add_proto_cluster(tmp_proto_clusters[proto_cluster_num])

            region.add_cand_cluster(cand_cluster)

        del tmp_proto_clusters, tmp_proto_cores, tmp_cand_clusters

    def __repr__(self) -> str:
        return f"GBK from {self.path} with {len(self.genes)} genes"
