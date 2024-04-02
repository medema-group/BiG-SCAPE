"""module containing integration tests for the distance generatio/comparison workflows"""

# from python
from unittest import TestCase
from pathlib import Path
from itertools import combinations

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation

# from other modules
import big_scape.data as bs_data
import big_scape.comparison as bs_comparison
import big_scape.file_input as bs_files
import big_scape.genbank as bs_gbk
from big_scape.genbank import (
    GBK,
    CDS,
    Region,
    CandidateCluster,
    ProtoCluster,
    ProtoCore,
)
import big_scape.enums as bs_enums


def create_mock_gbk(i, source_type: bs_enums.SOURCE_TYPE) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, 1, 0, 100, False, "test")
    gbk.region._db_id = i
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def create_mock_complete_single_gbk(
    i, source_type: bs_enums.SOURCE_TYPE, bgc_class
) -> GBK:
    region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
    region_feature.qualifiers = {
        "region_number": ["1"],
        "candidate_cluster_numbers": ["1"],
        "product": [bgc_class],
    }

    region = Region.parse_as5(region_feature)

    candidate_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
    candidate_cluster_feature.qualifiers = {
        "candidate_cluster_number": ["1"],
        "kind": ["single"],
        "protoclusters": ["1"],
        "product": [bgc_class],
    }

    candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

    protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature.qualifiers = {
        "protocluster_number": ["1"],
        "category": ["PKS"],
        "product": [bgc_class],
    }

    protocluster = ProtoCluster.parse(protocluster_feature)

    protocore_feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature.qualifiers = {
        "protocluster_number": ["1"],
        "product": [bgc_class],
    }

    protocore = ProtoCore.parse(protocore_feature)

    protocluster.add_proto_core(protocore)
    candidate_cluster.add_proto_cluster(protocluster)
    region.add_cand_cluster(candidate_cluster)

    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = region
    gbk.region._db_id = i
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def gen_mock_edge_list(
    edge_gbks: list[bs_gbk.GBK],
    edge_creation_function,
) -> list[
    tuple[int, int, float, float, float, float, int, bs_comparison.ComparableRegion]
]:
    edges = []
    for gbk_a, gbk_b in combinations(edge_gbks, 2):
        if gbk_a.region is None or gbk_b.region is None:
            continue
        if gbk_a.region._db_id is None or gbk_b.region._db_id is None:
            continue

        edges.append(edge_creation_function(gbk_a.region._db_id, gbk_b.region._db_id))

    return edges


def create_mock_edge_entirely_similar(a_id, b_id):
    return (
        a_id,
        b_id,
        0.0,
        1.0,
        1.0,
        1.0,
        1,
        bs_comparison.ComparableRegion(
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            False,
        ),
    )


def create_mock_edge_entirely_different(a_id, b_id):
    return (
        a_id,
        b_id,
        1.0,
        0.0,
        0.0,
        0.0,
        1,
        bs_comparison.ComparableRegion(
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            False,
        ),
    )


class TestComparison(TestCase):
    """Test class for hmm search/scan workflow"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_get_edge_param(self):
        """Tests whether the get_edge_param_id function correctly returns the correct params"""
        bs_data.DB.create_in_mem()

        run_1 = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }
        weights_1 = "mix"

        alignment_mode = run_1["alignment_mode"]

        edge_param_id = bs_comparison.edge_params_query(alignment_mode, weights_1)

        # the db is empty so we dont find an entry
        self.assertEqual(edge_param_id, None)

        edge_param_id = bs_comparison.edge_params_insert(alignment_mode, weights_1)

        edge_param_id = edge_param_id[0]

        # this is an empty db, so the first edge_param_id should be 1
        self.assertEqual(edge_param_id, 1)

        run_2 = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.GLOBAL,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }
        weights_2 = "classify"

        edge_param_id = bs_comparison.get_edge_param_id(run_2, weights_2)

        self.assertEqual(edge_param_id, 2)

        # the db now has two entries, and the function will still find the firs one
        # with the correct id
        edge_param_id = bs_comparison.get_edge_param_id(run_1, weights_1)

        self.assertEqual(edge_param_id, 1)

    def test_generate_edges(self):
        pass

    def test_calculate_scores_pair(self):
        pass

    def test_calculate_distances_mix_bin(self):
        """Tests the distance calculation workflow for mix mode"""

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
            "record_type": bs_enums.RECORD_TYPE.REGION,
        }
        weights = "mix"

        gbks = [
            create_mock_complete_single_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS")
            for i in range(3)
        ]

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        self.assertEqual(len(list_bgc_records), 3)

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        self.assertEqual(edge_param_id, 1)

        mix_bin = bs_comparison.generate_mix_bin(
            list_bgc_records, edge_param_id, run["record_type"]
        )

        expected_repr = "Bin 'mix': 3 pairs from 3 BGC records"

        self.assertEqual(str(mix_bin), expected_repr)

        expected_num_pairs = 3
        actual_num_pairs = mix_bin.num_pairs()

        self.assertEqual(actual_num_pairs, expected_num_pairs)

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        expected_num_pairs = 3
        actual_num_pairs = missing_edge_bin.num_pairs()

        self.assertEqual(actual_num_pairs, expected_num_pairs)

        edges = gen_mock_edge_list(gbks, create_mock_edge_entirely_similar)

        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)
        num_pairs = missing_edge_bin.num_pairs()

        self.assertEqual(num_pairs, 0)

    # def test_calculate_distances_legacy_classify(self):
    #     pass

    # def test_calculate_distances_classify(self):
    #     pass

    # def test_calculate_distances_query(self):
    #     pass
