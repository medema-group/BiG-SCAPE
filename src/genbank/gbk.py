"""Module containing code to load and store GBK files"""

# from python
import logging
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


class GBK:
    """
    Class to describe a given GBK file

    Attributes:
        path: Path
        metadata: Dict[str, str]
        region: Region
        nt_seq: SeqRecord.seq
        genes: list[CDS]
    """

    def __init__(self, path) -> None:
        self.path = path
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None
        self.nt_seq: SeqRecord.seq = None
        self.genes: List[Optional[CDS]] = []

    def save(self, commit=True):
        """Stores this GBK in the database

        Arguments:
            commit: commit immediately after executing the insert query"""
        gbk_table = DB.metadata.tables["gbk"]
        insert_query = (
            gbk_table.insert()
            .values(path=str(self.path), nt_seq=str(self.nt_seq))
            .compile()
        )

        DB.execute(insert_query)

        if commit:
            DB.commit()

    def save_all(self):
        """Stores this GBK and its children in the database. Does not commit immediately

        this function never commits"""
        self.save(False)
        self.region.save_all()

    @classmethod
    def parse(cls, path: Path):
        """Parses a GBK file and returns a GBK object with all necessary information"""
        gbk = cls(path)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))
        gbk.nt_seq = record.seq

        tmp_cand_clusters = {}
        tmp_proto_clusters = {}
        tmp_proto_cores = {}

        # go through features, load into tmp dicts indexed by feature number
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                if gbk.region is not None:
                    # this should not happen, but just in case
                    # since we only have one region on an object
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(feature)
                gbk.region = region

            if feature.type == "cand_cluster":
                cand_cluster = CandidateCluster.parse(feature)
                tmp_cand_clusters[cand_cluster.number] = cand_cluster

            if feature.type == "protocluster":
                proto_cluster = ProtoCluster.parse(feature)
                tmp_proto_clusters[proto_cluster.number] = proto_cluster

            if feature.type == "proto_core":
                proto_core = ProtoCore.parse(feature)
                tmp_proto_cores[proto_core.number] = proto_core

            if feature.type == "CDS":
                cds = CDS.parse(feature)
                gbk.genes.append(cds)

        # add features to parent objects
        for proto_cluster_num, proto_cluster in tmp_proto_clusters.items():
            proto_cluster.add_proto_core(tmp_proto_cores[proto_cluster_num])

        for cand_cluster in tmp_cand_clusters.values():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                cand_cluster.add_proto_cluster(tmp_proto_clusters[proto_cluster_num])

            region.add_cand_cluster(cand_cluster)

        del tmp_proto_clusters, tmp_proto_cores, tmp_cand_clusters

        return gbk
