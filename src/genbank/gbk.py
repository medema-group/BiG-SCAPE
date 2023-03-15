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
from src.genbank.protocluster import Protocluster
from src.genbank.proto_core import Protocore


class GBK:
    """
    Class to describe a given GBK file

    Attributes:
        path: Path
        metadata: Dict[str, str]
        region: Region
    """

    def __init__(self, path) -> None:
        self.path = path
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None

    @classmethod
    def parse(cls, path: Path):
        """Parses a GBK file and returns a GBK object with all necessary information"""
        gbk = cls(path)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))

        tmp_cand_clusters = {}
        tmp_protoclusters = {}
        tmp_proto_cores = {}

        # go through features, load into tmp dicts indexed by feature number
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                if gbk.region is not None:
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(feature)
                gbk.region = region

            elif feature.type == "cand_cluster":
                cand_cluster = CandidateCluster.parse(feature)
                tmp_cand_clusters[cand_cluster.number] = cand_cluster

            elif feature.type == "protocluster":
                protocluster = Protocluster.parse(feature)
                tmp_protoclusters[protocluster.number] = protocluster

            elif feature.type == "proto_core":
                proto_core = Protocore.parse(feature)
                tmp_proto_cores[proto_core.number] = proto_core

        # add features to parent objects
        for protocluster_nr, protocluster in tmp_protoclusters.items():
            protocluster.add_proto_core(tmp_proto_cores[protocluster_nr])

        for cand_cluster in tmp_cand_clusters.values():
            for protocluster_nr in cand_cluster.protoclusters.keys():
                cand_cluster.add_protocluster(tmp_protoclusters[protocluster_nr])

            region.add_cand_cluster(cand_cluster)

        del tmp_protoclusters, tmp_proto_cores, tmp_cand_clusters

        return gbk
