"""Module containing code to load and store GBK files"""

# from python
from __future__ import annotations
import logging

# from enum import Enum
from pathlib import Path
from typing import Dict, Optional
import hashlib


# from dependencies
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

# from other modules
from big_scape.errors import InvalidGBKError
from big_scape.data import DB
from big_scape.enums import SOURCE_TYPE
from big_scape.cli.config import BigscapeConfig

# from this module
from .region import Region
from .candidate_cluster import CandidateCluster
from .proto_cluster import ProtoCluster, MergedProtoCluster
from .proto_core import ProtoCore
from .cds import CDS


# class SOURCE_TYPE(Enum):
#     QUERY = "query"
#     MIBIG = "mibig"
#     REFERENCE = "reference"


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

    def __init__(self, path, hash, source_type) -> None:
        self.path: Path = path
        self.hash: str = hash
        self.metadata: Dict[str, str] = {}
        self.region: Optional[Region] = None
        self.nt_seq: SeqRecord.seq = None
        self.genes: list[CDS] = []
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
            hsp (hmmer.hsp): The HSP to be added to this CDS
            overlap_cutoff (float, optional): cutoff threshold for overlap. Defaults to
            0.1

        """
        # Filter out CDS for this gbk where CDS coordinates overlap by a certain amount.
        # The longest CDS will be chosen and other CDS will be discarded. This is done
        # in order to remove isoforms of the same gene in organisms that perform
        # alternative splicing

        # if no cds added yet, just add and continue
        if len(self.genes) == 0:
            self.genes.append(new_cds)
            return

        for cds_idx, old_cds in enumerate(self.genes):
            # ignore if these cds are on a different strand
            # BiG-SCAPE 1.0 does not perform this check
            # if not legacy_mode:
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
            self.genes.append(new_cds)
            return

        # if we got through all of that without exiting the function, we never replaced
        # a CDS so add a new one here
        self.genes.append(new_cds)

    def save(self, commit=True) -> None:
        """Stores this GBK in the database

        this returns the id of the GBK row

        Arguments:
            commit: commit immediately after executing the insert query"""

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        organism = self.metadata["organism"]
        taxonomy = self.metadata["taxonomy"]
        description = self.metadata["description"]

        gbk_table = DB.metadata.tables["gbk"]
        insert_query = (
            gbk_table.insert()
            .prefix_with("OR REPLACE")
            .values(
                path=str(self.path),
                hash=str(self.hash),
                nt_seq=str(self.nt_seq),
                organism=organism,
                taxonomy=taxonomy,
                description=description,
            )
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

        if return_row is None:
            raise RuntimeError("No return value from insert query")

        self._db_id = return_row[0]

        # only now that we have handled the return we can commit
        if commit:
            DB.commit()

    def save_all(self) -> None:
        """Stores this GBK and its children in the database. Does not commit immediately

        this function never commits"""
        self.save(False)

        if self.region is not None:
            self.region.save_all()

        for cds in self.genes:
            cds.save(False)

    @staticmethod
    def load_all() -> list[GBK]:
        """Load all GBK, CDS and BGCRecord objects from the database

        Returns:
            list[GBK]: _description_
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        gbk_table = DB.metadata.tables["gbk"]

        select_query = (
            gbk_table.select()
            .add_columns(
                gbk_table.c.id,
                gbk_table.c.hash,
                gbk_table.c.path,
                gbk_table.c.nt_seq,
                gbk_table.c.organism,
                gbk_table.c.taxonomy,
                gbk_table.c.description,
            )
            .compile()
        )

        cursor_result = DB.execute(select_query)

        gbk_dict = {}
        for result in cursor_result.all():
            new_gbk = GBK(Path(result.path), result.hash, "")
            new_gbk._db_id = result.id
            new_gbk.nt_seq = result.nt_seq
            new_gbk.metadata["organism"] = result.organism
            new_gbk.metadata["taxonomy"] = result.taxonomy
            new_gbk.metadata["description"] = result.description
            gbk_dict[result.id] = new_gbk

        # load GBK regions. This will also populate all record levels below region
        # e.g. candidate cluster, protocore if they exist

        Region.load_all(gbk_dict)

        CDS.load_all(gbk_dict)

        return list(gbk_dict.values())

    @staticmethod
    def load_many(input_gbks: list[GBK]) -> list[GBK]:
        """Load a list of GBK objects from the database

        Args:
            gbk_ids (list[int]): list of ids of gbk to load

        Returns:
            list[GBK]: loaded GBK objects
        """

        input_gbk_hashes = [gbk.hash for gbk in input_gbks]

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        gbk_table = DB.metadata.tables["gbk"]
        select_query = (
            gbk_table.select()
            .add_columns(
                gbk_table.c.id,
                gbk_table.c.hash,
                gbk_table.c.path,
                gbk_table.c.nt_seq,
                gbk_table.c.organism,
                gbk_table.c.taxonomy,
                gbk_table.c.description,
            )
            .where(gbk_table.c.hash.in_(input_gbk_hashes))
            .compile()
        )

        cursor_result = DB.execute(select_query)

        gbk_dict = {}
        for result in cursor_result.all():
            new_gbk = GBK(Path(result.path), result.hash, "")
            new_gbk._db_id = result.id
            new_gbk.nt_seq = result.nt_seq
            new_gbk.metadata["organism"] = result.organism
            new_gbk.metadata["taxonomy"] = result.taxonomy
            new_gbk.metadata["description"] = result.description
            gbk_dict[result.id] = new_gbk

        # load GBK regions. This will also populate all record levels below region
        # e.g. candidate cluster, protocore if they exist

        Region.load_all(gbk_dict)

        CDS.load_all(gbk_dict)

        return list(gbk_dict.values())

    @staticmethod
    def get_as_version(gbk_seq_record: SeqRecord) -> str:
        """Get AS version from GBK record

        Args:
            gbk_seq_record (SeqRecord): gbk seqrecord

        Returns:
            str: antismash version
        """

        try:
            as_version: str = gbk_seq_record.annotations["structured_comment"][
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
        run: dict,
        cds_overlap_cutoff: Optional[float] = None,
    ) -> GBK:
        """Parses a GBK file and returns a GBK object with all necessary information

        Args:
            path (Path): path to gbk file
            source_type (SOURCE_TYPE): A string describing the source of the GBK
            cds_overlap_cutoff (Optional[float]): cds region overlap cutoff to use.
            Defaults to None

        Returns:
            GBK: GBK object
        """

        # get unique content hash
        f = open(path, "r")
        data = f.read()
        f.close()
        data = data.encode("utf-8")  # type: ignore
        hash = hashlib.sha256(data).hexdigest()  # type: ignore

        gbk = cls(path, hash, source_type)

        # get record. should only ever be one for Antismash GBK
        record: SeqRecord = next(SeqIO.parse(path, "genbank"))
        gbk.nt_seq = record.seq

        gbk.metadata["description"] = record.description

        if "organism" in record.annotations:
            if record.annotations["organism"] == "":
                gbk.metadata["organism"] = "Unknown"
            else:
                gbk.metadata["organism"] = record.annotations["organism"]
        else:
            gbk.metadata["organism"] = "Unknown"

        if "taxonomy" in record.annotations:
            taxonomy = record.annotations["taxonomy"]
            if len(taxonomy) == 1:
                gbk.metadata["taxonomy"] = taxonomy[0]
            if len(taxonomy) > 1:
                taxonomy = ";".join(taxonomy)
                gbk.metadata["taxonomy"] = taxonomy
            else:
                gbk.metadata["taxonomy"] = "Unknown"
        else:
            gbk.metadata["taxonomy"] = "Unknown"

        as_version = GBK.get_as_version(record)
        # check if found version number is valid
        if not as_version[0].isnumeric():
            logging.error("%s: Invalid antiSMASH version in GBK header", gbk.path)
            raise InvalidGBKError()
        gbk.as_version = as_version

        if int(as_version[0]) < 6 and run["classify"] and run["legacy_weights"]:
            logging.error(
                "The parameters --classify and --legacy_weights are not compatible with "
                "gbks produced by antiSMASH versions under v6. You can either remove "
                "--legacy_weights or use --legacy_classify instead."
            )
            raise InvalidGBKError()

        if int(as_version[0]) >= 5:
            gbk.parse_as5up(record, cds_overlap_cutoff)
        if int(as_version[0]) == 4:
            # if regions are missing, gbk will be assumed to be as4
            # here we need to see if the user set --gbk-force to true

            gbk.parse_as4(record, cds_overlap_cutoff, run["force_gbk"])

        return gbk

    def parse_as4(
        self,
        record: SeqRecord,
        cds_overlap_cutoff: Optional[float] = None,
        force_gbk: bool = False,
    ) -> None:
        """Parses a GBK record of AS version 4 and returns a GBK object with all
        necessary information

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
                    logging.error(
                        "%s: GBK file contains more than one region", self.path
                    )
                    raise InvalidGBKError()

                region = Region.parse_as4(feature, parent_gbk=self)
                self.region = region

            if feature.type == "CDS":
                cds = CDS.parse(feature, parent_gbk=self)

                if cds is None:
                    continue

                cds.orf_num = orf_num
                orf_num += 1

                if cds_overlap_cutoff is not None:
                    self.add_cds_overlap_filter(cds, cds_overlap_cutoff)
                    continue

                self.genes.append(cds)

        if self.region is not None:
            return

        # If no cluster feature was found and force-gbk is false, GBK is invalid
        if not force_gbk:
            logging.error(
                "%s: GBK file does not contain an antiSMASH cluster feature", self.path
            )
            raise InvalidGBKError()

        # at this point we need to make a region object from the whole GBK
        logging.warning(
            "%s: GBK file does not contain an antiSMASH region feature. "
            "Using --force_gbk, assuming AS4",
            self.path,
        )

        self.region = Region.parse_full_region(record, parent_gbk=self)

    def parse_as5up(
        self,
        record: SeqRecord,
        cds_overlap_cutoff: Optional[float] = None,
    ) -> None:
        """Parses a GBK record of AS versions 5 and up and returns a GBK object with all
        necessary information

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
                    logging.error(
                        "%s: GBK file contains more than one region", self.path
                    )
                    raise InvalidGBKError()

                region = Region.parse_as5(feature, parent_gbk=self)
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

                if cds is None:
                    continue

                cds.orf_num = orf_num
                orf_num += 1

                if cds_overlap_cutoff is not None:
                    self.add_cds_overlap_filter(cds, cds_overlap_cutoff)
                    continue

                self.genes.append(cds)

        # check correct GBK format
        if self.region is None:
            logging.error(
                "%s: GBK file does not contain an antiSMASH region feature", self.path
            )
            raise InvalidGBKError()
        if not tmp_cand_clusters and not self.source_type == SOURCE_TYPE.MIBIG:
            logging.warning(
                "%s: GBK file does not contain any cand_cluster features", self.path
            )

        # add features to parent objects
        for proto_cluster_num, proto_cluster in tmp_proto_clusters.items():
            if proto_cluster_num not in tmp_proto_cores:
                logging.error(
                    "%s: protocluster %s has missing proto_core feature",
                    self.path,
                    proto_cluster_num,
                )
                raise InvalidGBKError()
            proto_cluster.add_proto_core(tmp_proto_cores[proto_cluster_num])

        updated_tmp_proto_clusters: dict[int, ProtoCluster] = {}

        self.collapse_hybrids_in_cand_clusters(
            tmp_cand_clusters, tmp_proto_clusters, updated_tmp_proto_clusters
        )

        # now we can add the protoclusters to the candidate clusters
        for number, cand_cluster in tmp_cand_clusters.items():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                cand_cluster.add_proto_cluster(
                    updated_tmp_proto_clusters[proto_cluster_num]
                )

            region.add_cand_cluster(cand_cluster)

        del (
            tmp_proto_clusters,
            tmp_proto_cores,
            tmp_cand_clusters,
            updated_tmp_proto_clusters,
        )

    def collapse_hybrids_in_cand_clusters(
        self, tmp_cand_clusters, tmp_proto_clusters, updated_tmp_proto_clusters
    ) -> None:
        """for chemical hybrid and interleaved candidate clusters, merge protoclusters into
        a single MergedProtoCluster, and update protocluster ids/numbers to reflect the new
        unique set of protoclusters

        Args:
            tmp_cand_clusters (dict[int, cand_clusters)
            tmp_proto_clusters (dict[int, ProtoCluster)
            updated_tmp_proto_clusters (dict[int, ProtoCluster)
        """

        # keep track of which protoclusters have been merged and old id: new id
        merged_protocluster_ids: dict[int, int] = {}
        merged_tmp_proto_clusters: dict[int, ProtoCluster] = {}

        # first we check which protoclusters need to be merged and do the merging
        for cand_cluster_num, cand_cluster in tmp_cand_clusters.items():
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                if proto_cluster_num not in tmp_proto_clusters:
                    logging.error(
                        "%s: cand_cluster %s has missing protocluster %s feature",
                        self.path,
                        cand_cluster_num,
                        proto_cluster_num,
                    )
                    raise InvalidGBKError()

            if cand_cluster.kind in BigscapeConfig.MERGED_CAND_CLUSTER_TYPE:
                # merge all the protoclusters in this candidate cluster
                protoclusters = [
                    tmp_proto_clusters[number]
                    for number in cand_cluster.proto_clusters.keys()
                ]
                merged_protocluster = MergedProtoCluster.merge(protoclusters)
                merged_tmp_proto_clusters[merged_protocluster.number] = (
                    merged_protocluster
                )

                # update the protocluster old:new ids for the merged protoclusters of this cand_cluster
                for proto_cluster_num in cand_cluster.proto_clusters.keys():
                    merged_protocluster_ids[proto_cluster_num] = (
                        merged_protocluster.number
                    )

        # now we build a new version of the tmp_proto_clusters dict that contains the merged protoclusters
        # as well as protoclusters which did not need merging, with updated unique IDs/numbers

        for cand_cluster_num, cand_cluster in tmp_cand_clusters.items():
            # for each cand_cluster, we also need to create an updated version of the proto_cluster dict
            # with the new IDs/numbers
            updated_proto_cluster_dict: dict[int, Optional[ProtoCluster]] = {}
            for proto_cluster_num in cand_cluster.proto_clusters.keys():
                if proto_cluster_num in merged_protocluster_ids.keys():
                    # this protocluster has been merged, so we need to add it to
                    # the dict with its new protocluster number
                    new_proto_cluster_num = merged_protocluster_ids[proto_cluster_num]
                    updated_tmp_proto_clusters[new_proto_cluster_num] = (
                        merged_tmp_proto_clusters[new_proto_cluster_num]
                    )
                    updated_proto_cluster_dict[new_proto_cluster_num] = None
                else:
                    # protoclusters which have not been merged are added to the dict as is
                    updated_tmp_proto_clusters[proto_cluster_num] = tmp_proto_clusters[
                        proto_cluster_num
                    ]
                    updated_proto_cluster_dict[proto_cluster_num] = None
            # we need to do this since the cand_cluster.add_proto_cluster method expects the dict to be there
            cand_cluster.proto_clusters = updated_proto_cluster_dict

        del (
            merged_protocluster_ids,
            merged_tmp_proto_clusters,
        )

    def __repr__(self) -> str:
        return f"GBK {self.path.name}, {len(self.genes)} genes"

    def __hash__(self) -> int:
        return hash(self.hash)

    def __eq__(self, other) -> bool:
        if not isinstance(other, GBK):
            return False

        if self.hash is None or other.hash is None:
            return False

        return self.hash == other.hash
