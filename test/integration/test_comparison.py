"""module containing integration tests for the distance generatio/comparison workflows"""

# from python
from unittest import TestCase
from pathlib import Path
from itertools import combinations

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation
from sqlalchemy import select

# from other modules
import big_scape.data as bs_data
import big_scape.comparison as bs_comparison
import big_scape.file_input as bs_files
import big_scape.genbank as bs_gbk
import big_scape.hmm as bs_hmm
from big_scape.genbank import (
    GBK,
    CDS,
    Region,
    CandidateCluster,
    ProtoCluster,
    ProtoCore,
)
import big_scape.distances.mix as bs_mix
import big_scape.distances.legacy_classify as bs_legacy_classify
import big_scape.distances.classify as bs_classify

import big_scape.enums as bs_enums


def create_mock_gbk(i, source_type: bs_enums.SOURCE_TYPE) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    add_mock_hsp_cds(cds)
    add_mock_hsp_alignment_hsp(cds.hsps[0])
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, 1, 0, 100, False, "test")
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def create_mock_complete_single_gbk(
    i, source_type: bs_enums.SOURCE_TYPE, bgc_class, bgc_category
) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)

    region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
    region_feature.qualifiers = {
        "region_number": ["1"],
        "candidate_cluster_numbers": ["1"],
        "product": [bgc_class],
    }

    region = Region.parse_as5(region_feature, parent_gbk=gbk)

    candidate_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
    candidate_cluster_feature.qualifiers = {
        "candidate_cluster_number": ["1"],
        "kind": ["single"],
        "protoclusters": ["1"],
        "product": [bgc_class],
    }

    candidate_cluster = CandidateCluster.parse(
        candidate_cluster_feature, parent_gbk=gbk
    )

    protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature.qualifiers = {
        "protocluster_number": ["1"],
        "category": [bgc_category],
        "product": [bgc_class],
    }

    protocluster = ProtoCluster.parse(protocluster_feature, parent_gbk=gbk)

    protocore_feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature.qualifiers = {
        "protocluster_number": ["1"],
        "product": [bgc_class],
    }

    protocore = ProtoCore.parse(protocore_feature, parent_gbk=gbk)

    protocluster.add_proto_core(protocore)
    candidate_cluster.add_proto_cluster(protocluster)
    region.add_cand_cluster(candidate_cluster)

    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = region
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def create_mock_complete_chemhyb_gbk(
    i, source_type: bs_enums.SOURCE_TYPE, bgc_class
) -> GBK:
    bgc_categories = bgc_class.split(".")

    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)

    region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
    region_feature.qualifiers = {
        "region_number": ["1"],
        "candidate_cluster_numbers": ["1"],
        "product": [bgc_class],
    }

    region = Region.parse_as5(region_feature, parent_gbk=gbk)

    candidate_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
    candidate_cluster_feature.qualifiers = {
        "candidate_cluster_number": ["1"],
        "kind": ["single"],
        "protoclusters": ["1", "2"],
        "product": [bgc_class],
    }

    candidate_cluster = CandidateCluster.parse(
        candidate_cluster_feature, parent_gbk=gbk
    )

    protocluster_feature_1 = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature_1.qualifiers = {
        "protocluster_number": ["1"],
        "category": [bgc_categories[0]],
        "product": [bgc_categories[0]],
    }

    protocluster_1 = ProtoCluster.parse(protocluster_feature_1, parent_gbk=gbk)

    protocore_feature_1 = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature_1.qualifiers = {
        "protocluster_number": ["1"],
        "product": [bgc_categories[0]],
    }

    protocore_1 = ProtoCore.parse(protocore_feature_1, parent_gbk=gbk)

    protocluster_feature_2 = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature_2.qualifiers = {
        "protocluster_number": ["2"],
        "category": [bgc_categories[1]],
        "product": [bgc_categories[1]],
    }

    protocluster_2 = ProtoCluster.parse(protocluster_feature_2, parent_gbk=gbk)

    protocore_feature_2 = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature_2.qualifiers = {
        "protocluster_number": ["2"],
        "product": [bgc_categories[1]],
    }

    protocore_2 = ProtoCore.parse(protocore_feature_2, parent_gbk=gbk)

    protocluster_1.add_proto_core(protocore_1)
    protocluster_2.add_proto_core(protocore_2)
    candidate_cluster.add_proto_cluster(protocluster_1)
    candidate_cluster.add_proto_cluster(protocluster_2)
    region.add_cand_cluster(candidate_cluster)

    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = region
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def add_mock_hsp_cds(cds: CDS, profile_accession: str = "PF01234") -> None:
    hsp = bs_hmm.HSP(cds, profile_accession, 1.0, 0, 100)
    cds.hsps.append(hsp)


def add_mock_hsp_alignment_hsp(hsp: bs_hmm.HSP, alignment: str = "AAAAAAAAAA") -> None:
    hsp_alignment = bs_hmm.HSPAlignment(hsp, alignment)
    hsp.alignment = hsp_alignment


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

    def test_fetch_records_from_db(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }
        weights = "mix"

        gbks_with_hsp = [
            create_mock_gbk(
                i,
                bs_enums.SOURCE_TYPE.QUERY,
            )
            for i in range(0, 3)
        ]

        gbks_no_hsp = [
            create_mock_complete_single_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "PKS", "PKS")
            for i in range(3, 6)
        ]

        gbks = gbks_with_hsp + gbks_no_hsp

        for gbk in gbks_with_hsp:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            gbk.genes[0].hsps[0].alignment.save()

        for gbk in gbks_no_hsp:
            gbk.save_all()

        bs_data.DB.commit()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        mix_bin = bs_comparison.generate_mix_bin(
            list_bgc_records, edge_param_id, run["record_type"]
        )

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        pair_ids = missing_edge_bin.generate_pair_ids()
        records = bs_comparison.workflow.fetch_records_from_database(pair_ids)

        self.assertEqual(len(records), 6)

    def test_calculate_scores_pair(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }
        weights = "mix"

        gbks_with_hsp = [
            create_mock_gbk(
                i,
                bs_enums.SOURCE_TYPE.QUERY,
            )
            for i in range(0, 3)
        ]

        gbks_no_hsp = [
            create_mock_complete_single_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "PKS", "PKS")
            for i in range(3, 6)
        ]

        gbks = gbks_with_hsp + gbks_no_hsp

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        mix_bin = bs_comparison.generate_mix_bin(
            list_bgc_records, edge_param_id, run["record_type"]
        )

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        pair_data = missing_edge_bin.generate_pairs()

        batch = bs_comparison.workflow.batch_generator(pair_data, 15)

        scores = bs_comparison.workflow.calculate_scores_pair(
            (
                batch,
                run["alignment_mode"],
                missing_edge_bin.edge_param_id,
                missing_edge_bin.weights,
            )
        )

        # 6 records make up a total of 15 edges
        self.assertEqual(len(scores), 15)

        distances_seen = []
        for score in scores:
            score = score[2]
            distances_seen.append(score)

        # 3 records are identical, so we expect 3 identical scores of 1.0
        # remaining are scores of 0.0 since the second batch of records do
        # not have any HSPs
        self.assertIn(1.0, distances_seen)
        self.assertIn(0.0, distances_seen)

    def test_generate_and_save_edges_workflow(self):
        """Tests the edge generation and saving workflow"""

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }
        weights = "mix"

        gbks = [
            create_mock_complete_single_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "PKS", "PKS")
            for i in range(3)
        ]

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        mix_bin = bs_comparison.generate_mix_bin(
            list_bgc_records, edge_param_id, run["record_type"]
        )

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        return_edges = []

        def callback(edges):
            nonlocal return_edges
            for edge in edges:
                return_edges.append(edge)

        bs_comparison.generate_edges(
            missing_edge_bin,
            run["alignment_mode"],
            run["cores"],
            run["cores"] * 2,
            callback,
            batch_size=5,
        )

        # check if the edges are correct
        # we expect 3 totally dissimilar edges,
        # as the regions are identical but have no HSPs

        self.assertEqual(len(return_edges), 3)

        scores = []
        for edge in return_edges:
            score = edge[2]
            score = round(score, 1)
            scores.append(score)

        self.assertEqual(scores, [1, 1, 1])

        # check that there are no edges in db
        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 0)

        bs_comparison.save_edges_to_db(return_edges, commit=True)

        # check that the edges are now in the db
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 3)

    def test_calculate_distances_mix_bin_worflow(self):
        """Tests the distance calculation workflow for mix mode"""

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
        }
        weights = "mix"

        gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
            )
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

    def test_calculate_distances_mix(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }

        gbks = [create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY) for i in range(3)]

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        # check that there are no edges in db
        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 0)

        bs_mix.calculate_distances_mix(run, list_bgc_records)

        # check that the edges are now in the db
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 3)

    def test_generate_bins_classify_worflow(self):
        """Tests the distance calculation workflow for classify mode"""

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": True,
        }

        pks_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
            )
            for i in range(0, 3)
        ]
        nrps_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS", "NRPS"
            )
            for i in range(3, 6)
        ]
        nrps_pks_gbks = [
            create_mock_complete_chemhyb_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "NRPS.PKS")
            for i in range(6, 9)
        ]

        gbks = pks_gbks + nrps_gbks + nrps_pks_gbks

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        self.assertEqual(len(list_bgc_records), 9)

        as_class_bins = list(
            bs_comparison.as_class_bin_generator(list_bgc_records, run)
        )

        self.assertEqual(len(as_class_bins), 3)

        bin_labels = [bin.label for bin in as_class_bins]

        self.assertEqual(bin_labels, ["T1PKS", "NRPS", "PKS"])

        bin_weights = [bin.weights for bin in as_class_bins]

        self.assertEqual(bin_weights, ["T1PKS", "NRPS", "PKSother"])

        bin_edge_param_ids = [bin.edge_param_id for bin in as_class_bins]

        self.assertEqual(bin_edge_param_ids, [1, 2, 3])

        bin_pairs = [bin.num_pairs() for bin in as_class_bins]

        self.assertEqual(bin_pairs, [3, 15, 3])

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": False,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": False,
        }

        as_class_bins = list(
            bs_comparison.as_class_bin_generator(list_bgc_records, run)
        )

        self.assertEqual(len(as_class_bins), 3)

        bin_labels = [bin.label for bin in as_class_bins]

        self.assertEqual(bin_labels, ["T1PKS", "NRPS", "NRPS.PKS"])

        bin_weights = [bin.weights for bin in as_class_bins]

        self.assertEqual(bin_weights, ["mix", "mix", "mix"])

        bin_edge_param_ids = [bin.edge_param_id for bin in as_class_bins]

        self.assertEqual(bin_edge_param_ids, [4, 4, 4])

        bin_pairs = [bin.num_pairs() for bin in as_class_bins]

        self.assertEqual(bin_pairs, [3, 3, 3])

    def test_generate_bins_legacy_classify_worflow(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": True,
        }

        pks_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
            )
            for i in range(0, 3)
        ]
        nrps_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS", "NRPS"
            )
            for i in range(3, 6)
        ]
        nrps_pks_gbks = [
            create_mock_complete_chemhyb_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS.otherks"
            )
            for i in range(6, 9)
        ]

        gbks = pks_gbks + nrps_gbks + nrps_pks_gbks

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        self.assertEqual(len(list_bgc_records), 9)

        as_class_bins = list(bs_comparison.legacy_bin_generator(list_bgc_records, run))

        self.assertEqual(len(as_class_bins), 3)

        bin_labels = [bin.label for bin in as_class_bins]

        self.assertEqual(bin_labels, ["PKSI", "NRPS", "PKSother"])

        bin_weights = [bin.weights for bin in as_class_bins]

        self.assertEqual(bin_weights, ["PKSI", "NRPS", "PKSother"])

        bin_edge_param_ids = [bin.edge_param_id for bin in as_class_bins]

        self.assertEqual(bin_edge_param_ids, [1, 2, 3])

        bin_pairs = [bin.num_pairs() for bin in as_class_bins]

        self.assertEqual(bin_pairs, [3, 15, 3])

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": False,
        }

        as_class_bins = list(bs_comparison.legacy_bin_generator(list_bgc_records, run))

        self.assertEqual(len(as_class_bins), 3)

        bin_labels = [bin.label for bin in as_class_bins]

        self.assertEqual(bin_labels, ["PKSI", "NRPS", "PKS-NRP_Hybrids"])

        bin_weights = [bin.weights for bin in as_class_bins]

        self.assertEqual(bin_weights, ["PKSI", "NRPS", "PKS-NRP_Hybrids"])

        bin_edge_param_ids = [bin.edge_param_id for bin in as_class_bins]

        self.assertEqual(bin_edge_param_ids, [1, 2, 4])

        bin_pairs = [bin.num_pairs() for bin in as_class_bins]

        self.assertEqual(bin_pairs, [3, 3, 3])

    def test_calculate_distances_legacy_classify(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": False,
            "cores": 1,
            "legacy_classify": True,
        }

        pks_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
            )
            for i in range(0, 3)
        ]
        nrps_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS", "NRPS"
            )
            for i in range(3, 6)
        ]
        nrps_pks_gbks = [
            create_mock_complete_chemhyb_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS.otherks"
            )
            for i in range(6, 9)
        ]

        gbks = pks_gbks + nrps_gbks + nrps_pks_gbks

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        self.assertEqual(len(list_bgc_records), 9)

        # bin 1 PKSI 3 pairs
        # bin 2 NRPS 3 pairs from 3 BGC records
        # bin PKS-NRP_Hybrids 3 pairs from 3 BGC records

        bs_legacy_classify.calculate_distances_legacy_classify(run, list_bgc_records)

        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 9)

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": True,
            "cores": 1,
            "legacy_classify": True,
        }

        # bin 1 PKSI 3 pairs -> those were already in db
        # bin 2 NRPS 15 pairs from 6 BGC records -> 3 pairs were already in db, so 12 new
        # bin PKSother 3 pairs from 3 BGC records -> 3 new pairs
        # 9 + 12 + 3 = 24

        bs_legacy_classify.calculate_distances_legacy_classify(run, list_bgc_records)

        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, 24)

    def test_calculate_distances_classify(self):
        """Tests the distance calculation workflow for classify mode"""

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": True,
            "cores": 1,
        }

        pks_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
            )
            for i in range(0, 3)
        ]
        nrps_gbks = [
            create_mock_complete_single_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS", "NRPS"
            )
            for i in range(3, 6)
        ]
        nrps_pks_gbks = [
            create_mock_complete_chemhyb_gbk(
                i, bs_enums.SOURCE_TYPE.QUERY, "NRPS.PKSother"
            )
            for i in range(6, 9)
        ]

        gbks = pks_gbks + nrps_gbks + nrps_pks_gbks

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        bs_classify.calculate_distances_classify(run, list_bgc_records)

        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        # bin 1 T1PKS has 3 pairs from 3 records
        # bin 2 NRPS has 15 pairs from 6 records
        # bin 3 PKSother has 3 pairs from 3 records

        self.assertEqual(rows, 21)

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": False,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": False,
            "cores": 1,
        }

        bs_classify.calculate_distances_classify(run, list_bgc_records)

        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        # bin 1 T1PKS has 3 pairs from 3 records with mix weights -> 3 new pairs
        # bin 2 NRPS has 3 pairs from 3 records with mix weights -> 3 new pairs
        # bin 3 NRPS.PKSother has 3 pairs from 3 records with mix weights -> 3 new pairs

        self.assertEqual(rows, 30)

    def test_test_generate_bins_query_workflow(self):
        self.skipTest("Not implemented")

    def test_calculate_distances_query(self):
        self.skipTest("Not implemented")
