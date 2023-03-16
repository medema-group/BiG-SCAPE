"""Module containing code to load and store GBK files"""

import logging
from pathlib import Path
from typing import Dict, Optional

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from src.errors.genbank import InvalidGBKError
from src.genbank.region import Region
from src.genbank.cand_cluster import CandidateCluster
from src.genbank.proto_cluster import Protocluster
from src.genbank.proto_core import Protocore


class GBK:
    """
    Class to describe a given GBK file

    Attributes:
        path: Path
        metadata: Dict[str, str]
        region: Region
        nt_seq: SeqRecord.seq
    """

    def __init__(self, path) -> None:
        self.path = path
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None
        self.nt_seq: SeqRecord.seq = None

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
                proto_cluster = Protocluster.parse(feature)
                tmp_proto_clusters[proto_cluster.number] = proto_cluster

            if feature.type == "proto_core":
                proto_core = Protocore.parse(feature)
                tmp_proto_cores[proto_core.number] = proto_core

        # add features to parent objects
        for proto_cluster_num, proto_cluster in tmp_proto_clusters.items():
            proto_cluster.add_proto_core(tmp_proto_cores[proto_cluster_num])

        for cand_cluster in tmp_cand_clusters.values():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                cand_cluster.add_proto_cluster(tmp_proto_clusters[proto_cluster_num])

            region.add_cand_cluster(cand_cluster)

        del tmp_proto_clusters, tmp_proto_cores, tmp_cand_clusters

        return gbk
