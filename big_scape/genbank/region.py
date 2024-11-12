"""Module containing code to load and store AntiSMASH regions"""

# from python
from __future__ import annotations
import logging
from typing import Dict, Optional, TYPE_CHECKING, Generator

# from dependencies
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

# from other modules
from big_scape.data import DB
from big_scape.errors import InvalidGBKError, InvalidGBKRegionChildError
from big_scape.enums import SOURCE_TYPE

# from this module
from .bgc_record import BGCRecord
from .candidate_cluster import CandidateCluster
from .cds import CDS


# from circular imports
if TYPE_CHECKING:  # pragma: no cover
    from big_scape.genbank import GBK  # imported earlier in file_input.load_files


class Region(BGCRecord):
    """
    Class to describe a region within an Antismash GBK

    Attributes:
        parent_gbk: GBK | None
        number: int
        contig_edge: Bool
        nt_start: int
        nt_stop: int
        product: str
        _db_id: int | None
        _families: dict[float, int]
        cand_clusters: Dict{number: int, CandidateCluster}
    """

    def __init__(
        self,
        parent_gbk: Optional[GBK],
        number: int,
        nt_start: int,
        nt_stop: int,
        contig_edge: Optional[bool],
        product: str,
    ):
        super().__init__(
            parent_gbk,
            number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )
        self.cand_clusters: Dict[int, Optional[CandidateCluster]] = {}

    def add_cand_cluster(self, cand_cluster: CandidateCluster) -> None:
        """Add a candidate cluster object to this region

        Args:
            cand_cluster (CandidateCluster): candidate cluster object

        Raises:
            InvalidGBKRegionChildError: Invalid gbk region child
        """

        if cand_cluster.number not in self.cand_clusters:
            raise InvalidGBKRegionChildError()

        self.cand_clusters[cand_cluster.number] = cand_cluster

    def save(self, commit=True) -> None:
        """Stores this region in the database

        Args:
            commit: commit immediately after executing the insert query"""
        super().save_record("region", commit)

    def save_all(self) -> None:
        """Stores this Region and its children in the database. Does not commit
        immediately
        """
        self.save(False)
        for candidate_cluster in self.cand_clusters.values():
            if candidate_cluster is None:
                continue

            if self._db_id is None:
                raise RuntimeError("Region has no database id!")

            candidate_cluster.save_all(self._db_id)

    @classmethod
    def parse_as5(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None) -> Region:
        """Creates a region object from a region feature in a GBK file

        Args:
            feature (SeqFeature): region(as5+) GBK feature

        Raises:
            InvalidGBKError: Invalid or missing fields

        Returns:
            Region: region object
        """
        err_path = parent_gbk.path if parent_gbk else ""

        if feature.type != "region":
            logging.error(
                "%s: feature is not of correct type! (expected: region, was: %s)",
                err_path,
                feature.type,
            )
            raise InvalidGBKError()

        # AS5 and up gbks have region features, as well as candidate clusters and
        # children classes (protocluster, protocore)

        if "region_number" not in feature.qualifiers:
            logging.error(
                "%s: region number qualifier not found in region feature!",
                err_path,
            )
            raise InvalidGBKError()

        region_number = int(feature.qualifiers["region_number"][0])

        nt_start, nt_stop, contig_edge = BGCRecord.parse_common(feature)

        # record may have multiple products. handle them here
        product = BGCRecord.parse_products(feature)

        # now we have all the parameters needed to assemble the region
        region = cls(
            parent_gbk,
            region_number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )

        if "candidate_cluster_numbers" not in feature.qualifiers:
            if parent_gbk is not None and parent_gbk.source_type == SOURCE_TYPE.MIBIG:
                # we know that MIBiG BGCs, although processed with versions of AS7 and
                # above dont always have features beyond region
                return region

            logging.error(
                "%s: candidate_cluster_numbers qualifier not found in region feature! "
                "Consider checking whether there is something special about this gbk",
                err_path,
            )
            raise InvalidGBKError()

        cand_clusters: dict[int, Optional[CandidateCluster]] = {}
        for cand_cluster_number in feature.qualifiers["candidate_cluster_numbers"]:
            cand_clusters[int(cand_cluster_number)] = None

        region.cand_clusters = cand_clusters
        return region

    @classmethod
    def parse_as4(cls, feature: SeqFeature, parent_gbk: Optional[GBK] = None) -> Region:
        """Creates a region object from a region feature in a GBK file

        Args:
            feature (SeqFeature): cluster (as4) GBK feature

        Raises:
            InvalidGBKError: Invalid or missing fields

        Returns:
            Region: region object
        """
        err_path = parent_gbk.path if parent_gbk else ""

        if feature.type != "cluster":
            logging.error(
                "%s: feature is not of correct type! (expected: cluster, was: %s)",
                err_path,
                feature.type,
            )
            raise InvalidGBKError()

        # AS4 gbks have cluster features instead of region, and no children features
        # we artifically input the info in the cluster feature into the Region object

        if (
            "note" not in feature.qualifiers
            or "Cluster number" not in feature.qualifiers["note"][0]
        ):
            logging.error(
                "%s: Cluster number qualifier not found in cluster feature!", err_path
            )
            raise InvalidGBKError()

        cluster_note_number = feature.qualifiers["note"][0]
        region_number = int(cluster_note_number.split(": ")[1])

        cand_clusters: dict[int, Optional[CandidateCluster]] = {}

        nt_start, nt_stop, contig_edge = BGCRecord.parse_common(feature)

        product = BGCRecord.parse_products_as4(feature)

        # now we have all the parameters needed to assemble the region
        region = cls(
            parent_gbk,
            region_number,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )
        region.cand_clusters = cand_clusters
        return region

    @classmethod
    def parse_full_region(cls, record: SeqRecord, parent_gbk: GBK) -> Region:
        """Generates a region from a full feature, without trying to parse any region feature

        Args:
            record (SeqRecord): SeqRecord object
            parent_gbk (GBK): parent GBK

        Returns:
            Region: region object
        """
        nt_start = 1
        nt_stop = len(record.seq) // 3
        contig_edge = False

        # record may have multiple products. handle them here
        product = "other"

        # now we have all the parameters needed to assemble the region
        region = cls(
            parent_gbk,
            1,
            nt_start,
            nt_stop,
            contig_edge,
            product,
        )

        region.cand_clusters = {}
        return region

    def get_cds_with_domains(
        self, return_all=True, reverse=False
    ) -> Generator[CDS, None, None]:
        return super().get_cds_with_domains(return_all, reverse)

    def __repr__(self):
        return f"{self.parent_gbk} Region {self.number} {self.nt_start}-{self.nt_stop} "

    @staticmethod
    def load_all(gbk_dict: dict[int, GBK]) -> None:
        """Load all Region objects from the database

        This function populates the region objects in the GBKs provided in the input
        gbk_dict

        Args:
            gbk_dict (dict[int, GBK]): Dictionary of GBK objects with database ids as
            keys. Used for reassembling the hierarchy
        """

        if not DB.metadata:
            raise RuntimeError("DB.metadata is None")

        record_table = DB.metadata.tables["bgc_record"]

        region_select_query = (
            record_table.select()
            .add_columns(
                record_table.c.id,
                record_table.c.gbk_id,
                record_table.c.parent_id,
                record_table.c.record_type,
                record_table.c.record_number,
                record_table.c.contig_edge,
                record_table.c.nt_start,
                record_table.c.nt_stop,
                record_table.c.product,
            )
            .where(record_table.c.record_type == "region")
            .where(record_table.c.gbk_id.in_(gbk_dict.keys()))
            .compile()
        )

        cursor_result = DB.execute(region_select_query)

        region_dict = {}

        for result in cursor_result.all():
            gbk = gbk_dict[result.gbk_id]

            new_region = Region(
                gbk,
                result.record_number,
                result.nt_start,
                result.nt_stop,
                result.contig_edge,
                result.product,
            )

            new_region._db_id = result.id

            # add to parent GBK
            gbk_dict[result.gbk_id].region = new_region

            # add to dictionary
            region_dict[result.id] = new_region

        CandidateCluster.load_all(region_dict)
