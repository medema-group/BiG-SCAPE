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
import big_scape.distances.query as bs_query
import big_scape.enums as bs_enums
import big_scape.network.network as bs_network


def create_mock_gbk_hsp(i, source_type: bs_enums.SOURCE_TYPE) -> GBK:
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


def create_mock_gbk(i, source_type: bs_enums.SOURCE_TYPE, product: str = "test") -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, 1, 0, 100, False, product)
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


def create_mock_query_dataset(run):
    """Creates a mock query dataset with 3 query records"""

    query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY, "NRPS")
    cds = query_gbk.genes[0]
    hsp_1 = bs_hmm.HSP(cds, "PF01234", 1.0, 0, 50)
    cds.hsps.append(hsp_1)
    hsp_1_alignment = bs_hmm.HSPAlignment(hsp_1, "AAAAAAAAAA")
    hsp_1.alignment = hsp_1_alignment

    gbk_1 = create_mock_gbk(1, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS")
    cds = gbk_1.genes[0]
    hsp_1 = bs_hmm.HSP(cds, "PF01234", 1.0, 0, 50)
    cds.hsps.append(hsp_1)
    hsp_1_alignment = bs_hmm.HSPAlignment(hsp_1, "AAAAAAAAAA")
    hsp_1.alignment = hsp_1_alignment

    hsp_2 = bs_hmm.HSP(cds, "PF05678", 1.0, 50, 100)
    cds.hsps.append(hsp_2)
    hsp_2_alignment = bs_hmm.HSPAlignment(hsp_2, "BBBBBBBBBB")
    hsp_2.alignment = hsp_2_alignment

    gbk_2 = create_mock_gbk(2, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS")
    cds = gbk_2.genes[0]
    hsp_2 = bs_hmm.HSP(cds, "PF05678", 1.0, 50, 100)
    cds.hsps.append(hsp_2)
    hsp_2_alignment = bs_hmm.HSPAlignment(hsp_2, "BBBBBBBBBB")
    hsp_2.alignment = hsp_2_alignment

    gbk_3 = create_mock_gbk(3, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS")
    cds = gbk_3.genes[0]
    hsp_3 = bs_hmm.HSP(cds, "PF09101", 1.0, 50, 100)
    cds.hsps.append(hsp_3)
    hsp_alignment = bs_hmm.HSPAlignment(hsp_3, "CCCCCCCCCC")
    hsp_3.alignment = hsp_alignment

    gbk_4 = create_mock_gbk(4, bs_enums.SOURCE_TYPE.REFERENCE, "T1PKS")

    gbks = [query_gbk, gbk_1, gbk_2, gbk_3, gbk_4]

    for gbk in gbks:
        gbk.save_all()
        for hsp in gbk.genes[0].hsps:
            hsp.save()
            hsp.alignment.save()

    query_gbk = gbks.pop(0)

    query_record = bs_files.get_all_bgc_records(run, [query_gbk])
    query_record = query_record[0]

    list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

    return query_record, list_bgc_records


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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }
        weights_1 = "mix"
        alignment_mode = run_1["alignment_mode"]
        extend_strat = run_1["extend_strategy"]

        edge_param_id = bs_comparison.edge_params_query(
            alignment_mode, weights_1, extend_strat
        )

        # the db is empty so we dont find an entry
        self.assertEqual(edge_param_id, None)

        edge_param_id = bs_comparison.edge_params_insert(
            alignment_mode, weights_1, extend_strat
        )

        edge_param_id = edge_param_id[0]

        # this is an empty db, so the first edge_param_id should be 1
        self.assertEqual(edge_param_id, 1)

        run_2 = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.GLOBAL,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.GREEDY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }

        gbks_with_hsp = [
            create_mock_gbk_hsp(
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

        mix_bin = bs_comparison.generate_mix_bin(list_bgc_records, run)

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        pair_ids = missing_edge_bin.generate_pair_ids()
        records = bs_comparison.workflow.fetch_records_from_database(pair_ids)

        self.assertEqual(len(records), 6)

    def test_calculate_scores_pair(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }

        gbks_with_hsp = [
            create_mock_gbk_hsp(
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

        mix_bin = bs_comparison.generate_mix_bin(list_bgc_records, run)

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        pair_data = missing_edge_bin.generate_pairs()

        batch = bs_comparison.workflow.batch_generator(pair_data, 15)

        scores = bs_comparison.workflow.calculate_scores_pair(
            (
                batch,
                run["alignment_mode"],
                run["extend_strategy"],
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }

        gbks = [
            create_mock_complete_single_gbk(i, bs_enums.SOURCE_TYPE.QUERY, "PKS", "PKS")
            for i in range(3)
        ]

        for gbk in gbks:
            gbk.save_all()

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        mix_bin = bs_comparison.generate_mix_bin(list_bgc_records, run)

        missing_edge_bin = bs_comparison.MissingRecordPairGenerator(mix_bin)

        return_edges = []

        def callback(edges):
            nonlocal return_edges
            for edge in edges:
                return_edges.append(edge)

        bs_comparison.generate_edges(
            missing_edge_bin,
            run["alignment_mode"],
            run["extend_strategy"],
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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

        mix_bin = bs_comparison.generate_mix_bin(list_bgc_records, run)

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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": True,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
        }

        gbks = [create_mock_gbk_hsp(i, bs_enums.SOURCE_TYPE.QUERY) for i in range(3)]

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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
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

    def test_missing_query_pair_generator_first_iteration(self):
        """Tests whether the QueryMissingRecordPairGenerator correctly generates a set of
        pairs in the first iteration of a specific network
        """

        bs_data.DB.create_in_mem()

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]

        query_pair_generator = bs_comparison.QueryRecordPairGenerator("mix", 1, "mix")
        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        # the state at the network at this point should be that all query to ref pairs
        # have been generated and their distances calculated. they are then added to
        # the network

        # the idea is that in the first iteration, since all query -> refs pairs
        # are in network, the missing pair generator should not generate any pairs
        # if we dont call the cycle.records() method
        query_pair_generator.add_records(source_records)

        missing_pair_generator = bs_comparison.QueryMissingRecordPairGenerator(
            query_pair_generator
        )

        # we will want to do this manually
        for idx, ref_gbk in enumerate(ref_gbks):
            # mypy please leave me alone
            if ref_gbk.region is None:
                continue
            if ref_gbk.region._db_id is None:
                continue
            if query_gbk.region is None:
                continue
            if query_gbk.region._db_id is None:
                continue

            # lets say two (0 and 1) of these distances are entirely similar
            if idx < 2:
                bs_comparison.save_edge_to_db(
                    (
                        query_gbk.region._db_id,
                        ref_gbk.region._db_id,
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
                )
            else:
                # the other two (2 and 3) are entirely distant
                bs_comparison.save_edge_to_db(
                    (
                        query_gbk.region._db_id,
                        ref_gbk.region._db_id,
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
                )

        # so now we have a network where the query is has distances of 0.0 to two of the
        # reference records, and 1.0 to the other two reference records

        # in our first iteration, we should see the following pairs:
        # (query, ref1)
        # (query, ref2)
        # (query, ref3)
        # (query, ref4)

        expected_pairs_query = set(
            [
                (query_gbk.region, ref_gbks[0].region),
                (query_gbk.region, ref_gbks[1].region),
                (query_gbk.region, ref_gbks[2].region),
                (query_gbk.region, ref_gbks[3].region),
            ]
        )

        actual_pairs_query = set(list(query_pair_generator.generate_pairs()))

        self.assertEqual(expected_pairs_query, actual_pairs_query)

        actual_num_pairs_query = query_pair_generator.num_pairs()

        self.assertEqual(4, actual_num_pairs_query)

        actual_pairs_missing = set(list(missing_pair_generator.generate_pairs()))

        self.assertEqual(set(), actual_pairs_missing)

        actual_num_pairs = missing_pair_generator.num_pairs()

        self.assertEqual(0, actual_num_pairs)

    def test_query_generators_workflow(self):
        """Tests whether the RefTorefPairGenerator correctly generates a set of
        pairs in the second iteration of a specific network
        """

        bs_data.DB.create_in_mem()

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]

        query_pair_generator = bs_comparison.QueryRecordPairGenerator("mix", 1, "mix")
        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        # the state at the network at this point should be that all query to ref pairs
        # have been generated and their distances calculated. they are then added to
        # the network

        # the idea is that in the first iteration, this should compare all connected
        # ref nodes to all singleton ref nodes

        query_pair_generator.add_records(source_records)

        missing_pair_generator = bs_comparison.QueryMissingRecordPairGenerator(
            query_pair_generator
        )

        # we will want to do this manually
        for idx, ref_gbk in enumerate(ref_gbks):
            # mypy please leave me alone
            if ref_gbk.region is None:
                continue
            if ref_gbk.region._db_id is None:
                continue
            if query_gbk.region is None:
                continue
            if query_gbk.region._db_id is None:
                continue

            # lets say two (0 and 1) of these distances are entirely similar
            if idx < 2:
                bs_comparison.save_edge_to_db(
                    (
                        query_gbk.region._db_id,
                        ref_gbk.region._db_id,
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
                )
            else:
                # the other two (2 and 3) are entirely distant
                bs_comparison.save_edge_to_db(
                    (
                        query_gbk.region._db_id,
                        ref_gbk.region._db_id,
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
                )

        # so now we have a network where the query is has distances of 0.0 to two of the
        # reference records, and 1.0 to the other two reference records

        # in our first iteration, we should see the following pairs:
        # (query, ref1)
        # (query, ref2)
        # (query, ref3)
        # (query, ref4)

        expected_pairs_query = set(
            [
                (query_gbk.region, ref_gbks[0].region),
                (query_gbk.region, ref_gbks[1].region),
                (query_gbk.region, ref_gbks[2].region),
                (query_gbk.region, ref_gbks[3].region),
            ]
        )

        actual_pairs_query = set(list(query_pair_generator.generate_pairs()))

        self.assertEqual(expected_pairs_query, actual_pairs_query)

        actual_pairs_missing = set(list(missing_pair_generator.generate_pairs()))

        self.assertEqual(set(), actual_pairs_missing)

        query_pair_generator.cycle_records(1.0)

        # passing distances will be query -> ref_1 and query -> ref_2

        expected_working_query = set([ref_gbks[1].region, ref_gbks[0].region])

        self.assertEqual(
            query_pair_generator.working_query_records, expected_working_query
        )

        expected_working_refs = set([ref_gbks[2].region, ref_gbks[3].region])

        self.assertEqual(
            query_pair_generator.working_ref_records, expected_working_refs
        )

        # # now we should see all the pairs between the nodes in expected working query and the
        # # working ref nodes, since none are in db yet.
        # # so we should see the following pairs:
        # # (ref1, ref0)
        # # (ref1, ref3)
        # # (ref2, ref0)
        # # (ref2, ref3)

        self.assertEqual(missing_pair_generator.num_pairs(), 4)

        # we will connect ref1 and ref3
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[1].region._db_id,
                ref_gbks[3].region._db_id,
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
        )

        # we will generate distant edges for the others
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[0].region._db_id,
                ref_gbks[1].region._db_id,
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
        )
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[0].region._db_id,
                ref_gbks[2].region._db_id,
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
        )
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[0].region._db_id,
                ref_gbks[3].region._db_id,
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
        )
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[1].region._db_id,
                ref_gbks[2].region._db_id,
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
        )
        bs_comparison.save_edge_to_db(
            (
                ref_gbks[2].region._db_id,
                ref_gbks[3].region._db_id,
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
        )

        # now the networks topology is this:
        # ref0    ref1
        #   \     /|
        #    \   / |
        #     \ /  |
        #     query|
        #          |
        #          |
        #          |
        # ref2    ref3
        # (There are more edges in the database, but these are all the ones with
        # distances < 1.0)

        # since we now put these edges in the databse, we should not see any missing pairs
        self.assertEqual(missing_pair_generator.num_pairs(), 0)

        query_pair_generator.cycle_records(1.0)

        # passing distances will be ref_1 -> ref_3

        expected_done = set([query_gbk.region, ref_gbks[0].region, ref_gbks[1].region])

        self.assertEqual(query_pair_generator.done_records, expected_done)

        expected_working_query = set([ref_gbks[3].region])

        self.assertEqual(
            query_pair_generator.working_query_records, expected_working_query
        )

        expected_working_refs = set([ref_gbks[2].region])

        self.assertEqual(
            query_pair_generator.working_ref_records, expected_working_refs
        )

        # one pair is going to be generated if we dont check the database
        self.assertEqual(query_pair_generator.num_pairs(), 1)

        # since the last missing pair is already in the database, we should not see any
        self.assertEqual(missing_pair_generator.num_pairs(), 0)

    def test_generate_bins_query_workflow(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": False,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "skip_propagation": False,
            "gcf_cutoffs": [0.1, 0.7],
        }

        query_record, list_bgc_records = create_mock_query_dataset(run)

        query_records = bs_query.get_query_records(run, list_bgc_records, query_record)

        self.assertEqual(len(query_records), 4)

        weights = "mix"

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)

        self.assertEqual(edge_param_id, 1)

        # final network should be:
        # query_gbk <-> gbk_1 <-> gbk_2  gbk_3

        # generate initial query -> ref pairs
        query_bin = bs_comparison.QueryRecordPairGenerator(
            "Query", edge_param_id, weights
        )
        query_bin.add_records(query_records)

        # fetch any existing distances from database
        missing_edge_bin = bs_comparison.QueryMissingRecordPairGenerator(query_bin)

        # there are no edges yet in db, and 1 query record and 3 ref records
        # so 3 edges need to be generated
        #     query_gbk <-> gbk_1
        #     query_gbk <-> gbk_2
        #     query_gbk <-> gbk_3
        self.assertEqual(missing_edge_bin.num_pairs(), 3)

        bs_query.calculate_distances(run, missing_edge_bin)

        # here the iterations will make the following edges:
        #     query_gbk <-> gbk_1   distance < 1.0
        #     query_gbk <-> gbk_2   distance = 1.0
        #     query_gbk <-> gbk_3   distance = 1.0
        #     gbk_1 <-> gbk_2       distance < 1.0
        #     gbk_1 <-> gbk_3       distance = 1.0
        #     gbk_2 <-> gbk_3       distance = 1.0

        # # all edges have been saved to db
        self.assertEqual(missing_edge_bin.num_pairs(), 0)

        distance_table = bs_data.DB.metadata.tables["distance"]
        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())
        self.assertEqual(rows, 6)

        # ref_to_ref_bin = bs_comparison.RefToRefRecordPairGenerator(
        #     "Ref_Ref", edge_param_id, weights
        # )
        # ref_to_ref_bin.add_records(query_records)

        # # now we have 1 connected ref (gbk_1) and 2 singleton refs, so 2 edges need to be generated
        # #      gbk_1 <-> gbk_2 and gbk_1 <-> gbk_3
        # # in a second iteration, we will have 1 connected ref and 1 singleton ref, so 1 edge
        # # more will be generated
        # #       gbk_2 <-> gbk_3
        # self.assertEqual(ref_to_ref_bin.num_pairs(), 2)

        # # generating distances -> Connected Reference to Singleton Reference
        # bs_query.calculate_distances(
        #     run, ref_to_ref_bin, "Singleton Reference to Connected Reference"
        # )

        # self.assertEqual(ref_to_ref_bin.num_pairs(), 0)

        # select_statement = select(distance_table.c.distance)
        # rows = len(bs_data.DB.execute(select_statement).fetchall())
        # self.assertEqual(rows, 6)

        # now we make any last connected ref <-> connected ref pairs that are missing
        # get all the edges in the query connected component
        query_connected_component = bs_network.get_connected_components(
            1, edge_param_id, query_bin, 1
        )

        # get_connected_components returns a list of connected components, we only want the first one
        query_connected_component = next(query_connected_component)

        # there are only 2 edges that have scores < 1.0
        #   gbk_1 <-> gbk_2
        #   query_gbks <-> gbk_1
        self.assertEqual(len(query_connected_component), 2)

        query_nodes = bs_network.get_nodes_from_cc(
            query_connected_component, query_records
        )

        self.assertEqual(len(query_nodes), 3)

        self.assertEqual(query_record, query_nodes[0])
        self.assertEqual(query_records[1], query_nodes[1])
        self.assertEqual(query_records[2], query_nodes[2])

        query_connected_bin = bs_comparison.RecordPairGenerator(
            "Query", edge_param_id, record_type=run["record_type"]
        )
        query_connected_bin.add_records(query_nodes)

        # fetch any existing distances from database
        missing_ref_edge_bin = bs_comparison.MissingRecordPairGenerator(
            query_connected_bin
        )

        # there is 1 edge that is not in the cc if considering the cutoff
        #   gbk_1 <-> gbk_2
        # but this edge is already in database with a distance of 1!
        # (from the first iteration)
        self.assertEqual(missing_ref_edge_bin.num_pairs(), 0)

        # calculate distances -> Connected Reference to Connected Reference
        bs_query.calculate_distances(
            run,
            missing_ref_edge_bin,
        )

        select_statement = select(distance_table.c.distance)
        rows = len(bs_data.DB.execute(select_statement).fetchall())
        self.assertEqual(rows, 6)

    def test_calculate_distances_query(self):
        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": False,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "skip_propagation": False,
            "run_id": 1,
            "gcf_cutoffs": [0.1, 0.8],
        }

        query_record, list_bgc_records = create_mock_query_dataset(run)

        query_to_ref_bin = bs_comparison.QueryRecordPairGenerator("Query_Ref", 1, "mix")
        query_records = bs_query.get_query_records(run, list_bgc_records, query_record)

        query_to_ref_bin.add_records(query_records)

        bs_query.calculate_distances_query(run, list_bgc_records, query_record)

        query_connected_component = bs_network.get_connected_components(
            1, 1, query_to_ref_bin, 1
        )

        query_connected_component = next(query_connected_component)

        # there are only 2 edges that have scores < 1.0
        #   gbk_1 <-> gbk_2
        #   query_gbks <-> gbk_1
        self.assertEqual(len(query_connected_component), 2)

        query_records = [query_record] + list_bgc_records

        query_nodes = bs_network.get_nodes_from_cc(
            query_connected_component, query_records
        )

        self.assertEqual(len(query_nodes), 3)

        self.assertEqual(query_record, query_nodes[0])
        self.assertEqual(query_records[1], query_nodes[1])
        self.assertEqual(query_records[2], query_nodes[2])

    def test_get_query_records(self):
        bs_data.DB.create_in_mem()

        gbk_1 = create_mock_complete_single_gbk(
            1, bs_enums.SOURCE_TYPE.REFERENCE, "T1PKS", "PKS"
        )

        gbk_2 = create_mock_complete_single_gbk(
            2, bs_enums.SOURCE_TYPE.REFERENCE, "T2PKS", "PKS"
        )

        gbk_3 = create_mock_complete_single_gbk(
            3, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS", "NRPS"
        )

        gbk_4 = create_mock_complete_single_gbk(
            4, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS.T1PKS", "NRPS.PKS"
        )

        gbk_5 = create_mock_complete_single_gbk(
            5, bs_enums.SOURCE_TYPE.REFERENCE, "NRPS.T2PKS", "NRPS.PKS"
        )

        gbks = [gbk_1, gbk_2, gbk_3, gbk_4, gbk_5]

        for gbk in gbks:
            gbk.save_all()

        query_gbk = create_mock_complete_single_gbk(
            0, bs_enums.SOURCE_TYPE.QUERY, "T1PKS", "PKS"
        )
        query_gbk.save_all()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": False,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
        }

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        query_record = bs_files.get_all_bgc_records(run, [query_gbk])
        query_record = query_record[0]

        query_records = bs_query.get_query_records(run, list_bgc_records, query_record)

        self.assertEqual(len(query_records), 2)

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "extend_strategy": bs_enums.EXTEND_STRATEGY.LEGACY,
            "legacy_weights": False,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "cores": 1,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }

        list_bgc_records = bs_files.get_all_bgc_records(run, gbks)

        query_record = bs_files.get_all_bgc_records(run, [query_gbk])
        query_record = query_record[0]

        query_records = bs_query.get_query_records(run, list_bgc_records, query_record)

        self.assertEqual(len(query_records), 3)
