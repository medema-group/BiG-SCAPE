"""Module containing code to load and store GBK files"""

# from python
import logging
from enum import Enum
from pathlib import Path
from typing import Dict, Optional


# from dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
from sortedcontainers import SortedList

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
        genes: SortedList[CDS]
        as_version: str
        source_type: SOURCE_TYPE
    """

    def __init__(self, path, source_type) -> None:
        self.path: Path = path
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None
        self.nt_seq: SeqRecord.seq = None
        self.genes: SortedList[CDS] = SortedList()
        self.as_version: Optional[str] = None
        self.source_type: SOURCE_TYPE = source_type

        # db specific fields
        self._db_id: Optional[int] = None

    def add_cds_overlap_filter(self, new_cds: CDS, cds_overlap_cutoff=0.1) -> None:
        """Adds a CDS to this GBK. Performs overlap cutoff filtering by calculating the
        percentage overlap of the incoming CDS with other CDS in this GBK.

        If the percentage overlap is greater than the cutoff, this keeps whichever CDS
        is longer. If scores are equal, this keeps the HSP with the earliest
        start position. Bitscores are rounded to 1 decimal position when compared.

        The above behavior should mirror BiG-SCAPE 1.0 behavior

        Args:
            hsp (src.hmmer.hsp): The HSP to be added to this CDS
            overlap_cutoff (float, optional): cutoff threshold for overlap. Defaults to
            0.1

        """
        # Filter out CDS for this gbk where CDS coordinates overlap by a certain amount.
        # The longest CDS will be chosen and other CDS will be discarded. This is done
        # in order to remove isoforms of the same gene in organisms that perform
        # alternative splicing

        # if no cds added yet, just add and continue
        if len(self.genes) == 0:
            self.genes.add(new_cds)
            return

        for cds_idx, old_cds in enumerate(self.genes):
            # ignore if these cds are on a different strand
            if old_cds.strand != new_cds.strand:
                continue

            overlap_nt = CDS.len_nt_overlap(old_cds, new_cds)

            # ignore these cds if there is no overlap at all
            if overlap_nt == 0:
                continue

            overlap_aa = overlap_nt / 3

            # we determine overlap as a percentage of shortest cds, so get that
            old_cds_aa_len = len(old_cds.aa_seq)
            new_cds_aa_len = len(new_cds.aa_seq)
            shortest_aa_len = min(old_cds_aa_len, new_cds_aa_len)

            logging.debug(
                "%s, %d, %d-%d : %d-%d",
                self.path,
                overlap_nt,
                old_cds.nt_start,
                old_cds.nt_stop,
                new_cds.nt_start,
                new_cds.nt_stop,
            )

            # not over cutoff? skip
            if overlap_aa <= cds_overlap_cutoff * shortest_aa_len:
                continue

            # new cds is shorter? skip
            if new_cds_aa_len < old_cds_aa_len:
                logging.debug(
                    "Removing %s because it overlaps with another CDS", new_cds
                )
                return

            # overwrite if equal or larger, same as BiG-SCAPE 1.0
            logging.debug(
                "Removing %s because it overlaps with another CDS", self.genes[cds_idx]
            )
            del self.genes[cds_idx]
            self.genes.add(new_cds)
            return

        # if we got through all of that without exiting the function, we never replaced
        # a CDS so add a new one here
        self.genes.add(new_cds)

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

        for cds in self.genes:
            cds.save(False)

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
    def parse(
        cls,
        path: Path,
        source_type: SOURCE_TYPE,
        cds_overlap_cutoff: Optional[float] = None,
    ):
        """Parses a GBK file and returns a GBK object with all necessary information

        Args:
            path (Path): path to gbk file
            source_type (SOURCE_TYPE): A string describing the source of the GBK
            cds_overlap_cutoff (Optional[float]): cds region overlap cutoff to use.
            Defaults to None

        Returns:
            GBK: GBK object
        """

        gbk = cls(path, source_type)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))
        gbk.nt_seq = record.seq

        if "organism" in record.annotations:
            gbk.metadata["organism"] = record.annotations["organism"]

        as_version = GBK.get_as_version(record)
        gbk.as_version = as_version

        if int(as_version[0]) >= 5:
            gbk.parse_as5up(record, cds_overlap_cutoff)
        if int(as_version[0]) == 4:
            gbk.parse_as4(record, cds_overlap_cutoff)

        return gbk

    def parse_as4(self, record: SeqRecord, cds_overlap_cutoff: Optional[float] = None):
        """Parses a GBK record of AS version 4 and returns a GBK object with all necessary information

        Args:
            record (SeqRecord): gbk file record
            cds_overlap_cutoff (Optional[float]): cds region overlap cutoff to use.
            Defaults to None

        Raises:
            InvalidGBKError: Invalid or missing fields in gbk record
        """

        # go through features, load into tmp dicts indexed by feature number
        orf_num = 1
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "cluster":
                if self.region is not None:
                    # this should not happen, but just in case
                    # since we only have one region on an object
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(feature, parent_gbk=self)
                self.region = region

            if feature.type == "CDS":
                cds = CDS.parse(feature, parent_gbk=self)

                cds.orf_num = orf_num
                orf_num += 1

                if cds_overlap_cutoff is not None:
                    self.add_cds_overlap_filter(cds, cds_overlap_cutoff)
                    continue

                self.genes.add(cds)

    def parse_as5up(
        self, record: SeqRecord, cds_overlap_cutoff: Optional[float] = None
    ):
        """Parses a GBK record of AS versions 5 and up and returns a GBK object with all necessary information

        Args:
            record (SeqRecord): gbk file record
            cds_overlap_cutoff (Optional[float]): cds region overlap cutoff to use.
            Defaults to None

        Raises:
            InvalidGBKError: Invalid or missing fields in gbk record
        """

        tmp_cand_clusters = {}
        tmp_proto_clusters = {}
        tmp_proto_cores = {}

        # go through features, load into tmp dicts indexed by feature number
        orf_num = 1
        feature: SeqFeature
        for feature in record.features:
            if feature.type == "region":
                if self.region is not None:
                    # this should not happen, but just in case
                    # since we only have one region on an object
                    logging.error("GBK file provided contains more than one region")
                    raise InvalidGBKError()

                region = Region.parse(feature, parent_gbk=self)
                self.region = region

            if feature.type == "cand_cluster":
                cand_cluster = CandidateCluster.parse(feature, parent_gbk=self)
                tmp_cand_clusters[cand_cluster.number] = cand_cluster

            if feature.type == "protocluster":
                proto_cluster = ProtoCluster.parse(feature, parent_gbk=self)
                tmp_proto_clusters[proto_cluster.number] = proto_cluster

            if feature.type == "proto_core":
                proto_core = ProtoCore.parse(feature, parent_gbk=self)
                tmp_proto_cores[proto_core.number] = proto_core

            if feature.type == "CDS":
                cds = CDS.parse(feature, parent_gbk=self)

                cds.orf_num = orf_num
                orf_num += 1

                if cds_overlap_cutoff is not None:
                    self.add_cds_overlap_filter(cds, cds_overlap_cutoff)
                    continue

                self.genes.add(cds)

        # add features to parent objects
        for proto_cluster_num, proto_cluster in tmp_proto_clusters.items():
            proto_cluster.add_proto_core(tmp_proto_cores[proto_cluster_num])

        for cand_cluster in tmp_cand_clusters.values():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                cand_cluster.add_proto_cluster(tmp_proto_clusters[proto_cluster_num])

            region.add_cand_cluster(cand_cluster)

        del tmp_proto_clusters, tmp_proto_cores, tmp_cand_clusters

    def __repr__(self) -> str:
        return f"GBK {self.path.name}, {len(self.genes)} genes"

    def __hash__(self) -> int:
        return hash(self.path)
