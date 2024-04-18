"""Contains tests for the pair class and its methods"""

# from python
from unittest import TestCase
from pathlib import Path
from collections import OrderedDict

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation
from big_scape.comparison.record_pair import RecordPair

# from other modules
from big_scape.genbank import GBK, BGCRecord, CDS, Region, ProtoCluster, ProtoCore
from big_scape.comparison import (
    RecordPairGenerator,
    QueryToRefRecordPairGenerator,
    RefToRefRecordPairGenerator,
    ConnectedComponentPairGenerator,
    save_edge_to_db,
    get_record_category,
    get_legacy_weights_from_category,
    as_class_bin_generator,
)
from big_scape.comparison import generate_mix_bin
from big_scape.genbank.candidate_cluster import CandidateCluster

import big_scape.hmm as bs_hmm
import big_scape.data as bs_data
import big_scape.enums as bs_enums
import big_scape.comparison as bs_comparison


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


def add_mock_hsp_cds(cds: CDS) -> None:
    hsp = bs_hmm.HSP(cds, "PF01234", 1.0, 0, 100)
    cds.hsps.append(hsp)


def add_mock_hsp_alignment_hsp(hsp: bs_hmm.HSP) -> None:
    hsp_alignment = bs_hmm.HSPAlignment(hsp, "AAAAAAAAAA")
    hsp.alignment = hsp_alignment


class TestBGCPair(TestCase):
    """Contains tests to test pair objects used in binning"""

    def test_pair_repr(self):
        """Tests whether calling str() on a bin object returns an expected string
        representation of the object
        """
        gbk = GBK(Path("test"), "test", bs_enums.SOURCE_TYPE.QUERY)

        bgc_a = BGCRecord(gbk, 0, 0, 10, False, "")
        bgc_b = BGCRecord(gbk, 0, 10, 20, False, "")

        pair = RecordPair(bgc_a, bgc_b)

        expected_repr = (
            "Pair GBK test, 0 genes Record (superclass) 0-10 - GBK test, 0 genes "
            "Record (superclass) 10-20"
        )

        actual_repr = str(pair)

        self.assertEqual(expected_repr, actual_repr)

    def test_pair_no_parent_gbk(self):
        """Tests whether initialization of a BGC pair where one of the BGCs does not
        have a parent GBK correctly throws a ValueError
        """
        gbk = GBK("", "", "test")

        bgc_a = BGCRecord(gbk, 0, 0, 10, False, "")

        # b is missing GBK
        bgc_b = BGCRecord(None, 0, 0, 10, False, "")

        self.assertRaises(ValueError, RecordPair, bgc_a, bgc_b)


class TestBGCBin(TestCase):
    def clean_db(self):
        """Closes the database connection and resets the metadata"""
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_bin_repr(self):
        """Tests whether calling str() on a bin object returns an expected string
        representation of the object
        """
        parent_gbk = GBK(Path("test"), "test", source_type=bs_enums.SOURCE_TYPE.QUERY)
        bgc_a = BGCRecord(parent_gbk, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk, 0, 0, 10, False, "")

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = RecordPairGenerator("test", 1)

        new_bin.add_records(bgc_list)

        # expected representation of the bin object
        expected_repr = "Bin 'test': 3 pairs from 3 BGC records"

        actual_repr = str(new_bin)

        self.assertEqual(expected_repr, actual_repr)

    def test_num_pairs_too_few_records(self):
        """tests if bin.num_pairs() correctly returns 0 if there is only one record in the bin"""

        gbk_a = GBK(Path("test1.gbk"), "test", "test")
        bgc_a = BGCRecord(gbk_a, 0, 0, 10, False, "")

        new_bin = RecordPairGenerator("test", 1)

        new_bin.add_records([bgc_a])

        expected_num_pairs = 0
        actual_num_pairs = new_bin.num_pairs()

        self.assertEqual(expected_num_pairs, actual_num_pairs)

    def test_num_pairs_correct_with_query_ref(self):
        """Tests whether bin.num_pairs() correctly returns all query and ref but not ref <-> ref pairs"""

        parent_gbk_query = GBK(
            Path("test"), "test", source_type=bs_enums.SOURCE_TYPE.QUERY
        )
        parent_gbk_ref = GBK(
            Path("test"), "test", source_type=bs_enums.SOURCE_TYPE.REFERENCE
        )
        bgc_a = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_b = BGCRecord(parent_gbk_query, 0, 0, 10, False, "")
        bgc_c = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")
        bgc_d = BGCRecord(parent_gbk_ref, 0, 0, 10, False, "")

        bgc_list = [bgc_a, bgc_b, bgc_c, bgc_d]

        new_bin = QueryToRefRecordPairGenerator("test", 1, "mix")

        new_bin.add_records(bgc_list)

        expected_num_pairs = 5
        actual_num_pairs = new_bin.num_pairs()

        self.assertEqual(expected_num_pairs, actual_num_pairs)

    def test_num_pairs_corrected_multiple_subrecords(self):
        """Tests whether number of expected pairs is corrected for subrecords from the
        same gbk"""
        bs_data.DB.create_in_mem()

        gbk1 = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        pclust1 = ProtoCluster(gbk1, 1, 10, 100, False, "", {})
        pclust2 = ProtoCluster(gbk1, 2, 50, 200, False, "", {})
        gbk1.region.cand_clusters[1] = CandidateCluster(
            gbk1, 0, 10, 100, False, "", "", {1: pclust1, 2: pclust2}
        )
        gbk1.save_all()

        gbk2 = create_mock_gbk(1, bs_enums.SOURCE_TYPE.QUERY)
        pclust3 = ProtoCluster(gbk2, 1, 10, 100, False, "", {})
        pclust4 = ProtoCluster(gbk2, 2, 50, 200, False, "", {})
        pclust5 = ProtoCluster(gbk2, 3, 200, 300, False, "", {})
        gbk2.region.cand_clusters[1] = CandidateCluster(
            gbk1, 0, 10, 100, False, "", "", {1: pclust3, 2: pclust4, 3: pclust5}
        )
        gbk2.save_all()

        gbk3 = create_mock_gbk(2, bs_enums.SOURCE_TYPE.QUERY)
        gbk3.save_all()

        records = [pclust1, pclust2, pclust3, pclust4, pclust5, gbk3.region]
        bin = RecordPairGenerator(
            "test", 1, record_type=bs_enums.RECORD_TYPE.PROTO_CLUSTER
        )
        bin.add_records(records)

        # with 6 total records, we would normally expect 15 pairs
        # subtract one for gbk1 and subtract 3 for gbk2
        expected_pairs = 11
        actual_pairs = bin.num_pairs()

        self.assertEqual(expected_pairs, actual_pairs)

    def test_legacy_sorting(self):
        """Tests whether the legacy sorting option in bin.pairs() correctly orders the pairs"""

        gbk_a = GBK(Path("test1.gbk"), "test1", "test")
        bgc_a = BGCRecord(gbk_a, 0, 0, 10, False, "")
        bgc_a._db_id = 1
        gbk_b = GBK(Path("test2.gbk"), "test2", "test")
        bgc_b = BGCRecord(gbk_b, 0, 0, 10, False, "")
        bgc_b._db_id = 2
        gbk_c = GBK(Path("test3.gbk"), "test3", "test")
        bgc_c = BGCRecord(gbk_c, 0, 0, 10, False, "")
        bgc_c._db_id = 3

        # due to the order, this should generate a list of pairs as follows without legacy sort:
        # bgc_a, bgc_c
        # bgc_a, bgc_b
        # bgc_c, bgc_b
        bgc_list = [bgc_a, bgc_c, bgc_b]

        new_bin = RecordPairGenerator("test", 1)

        new_bin.add_records(bgc_list)

        # expected list should correctly sort the third entry int the list to be bgc_b, bgc_c
        expected_pair_list = [
            (bgc_a, bgc_c),
            (bgc_a, bgc_b),
            (bgc_b, bgc_c),
        ]

        actual_pair_list = [
            pair for pair in new_bin.generate_pairs(legacy_sorting=True)
        ]

        self.assertEqual(expected_pair_list, actual_pair_list)

    def test_generate_pairs(self):
        """Tests whether bin.generate_pairs() correctly generates all pairs"""

        gbk_a = GBK(Path("test1.gbk"), "test1", "test")
        bgc_a = BGCRecord(gbk_a, 0, 0, 10, False, "")
        bgc_a._db_id = 1
        gbk_b = GBK(Path("test2.gbk"), "test2", "test")
        bgc_b = BGCRecord(gbk_b, 0, 0, 10, False, "")
        bgc_b._db_id = 2
        gbk_c = GBK(Path("test3.gbk"), "test3", "test")
        bgc_c = BGCRecord(gbk_c, 0, 0, 10, False, "")
        bgc_c._db_id = 3

        # due to the order, this should generate a list of pairs as follows without legacy sort:
        # bgc_a, bgc_c
        # bgc_a, bgc_b
        # bgc_c, bgc_b
        bgc_list = [bgc_a, bgc_c, bgc_b]

        new_bin = RecordPairGenerator("test", 1)

        new_bin.add_records(bgc_list)

        actual_pair_list = [pair for pair in new_bin.generate_pairs()]

        expected_pair_ids = [
            (pair[0]._db_id, pair[1]._db_id) for pair in actual_pair_list
        ]

        actual_pair_ids = [pair for pair in new_bin.generate_pair_ids()]

        self.assertEqual(expected_pair_ids, actual_pair_ids)

    def test_query_to_ref_pair_generator(self):
        """Tests whether the QueryToRefPairGenerator correctly generates a set of
        pairs in a specific network
        """

        bs_data.DB.create_in_mem()

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 4)
        ]

        query_to_ref_pair_generator = QueryToRefRecordPairGenerator("mix", 1, "mix")
        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)

        query_to_ref_pair_generator.add_records(source_records)

        # expected edges    def test_ref_to_ref_pair_generator_first_iteration(self):

        expected_pairs = []
        for ref_gbk in ref_gbks:
            expected_pair = (query_gbk.region, ref_gbk.region)
            expected_pairs.append(expected_pair)

        # get all edges
        actual_pairs = list(query_to_ref_pair_generator.generate_pairs())

        self.assertListEqual(expected_pairs, actual_pairs)

    def test_ref_to_ref_pair_generator_first_iteration(self):
        """Tests whether the RefTorefPairGenerator correctly generates a set of
        pairs in the first iteration of a specific network
        """
        bs_data.DB.create_in_mem()

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]

        ref_to_ref_pair_generator = RefToRefRecordPairGenerator("mix", 1, "mix")
        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        # the state at the network at this point should be that all query to ref pairs
        # have been generated and their distances calculated. they are then added to
        # the network

        # the idea is that in the first iteration, this should compare all connected
        # ref nodes to all singleton ref nodes
        ref_to_ref_pair_generator.add_records(source_records)

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
                save_edge_to_db(
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
                save_edge_to_db(
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

        # so now we have a network where the query is connected to 2 of the reference
        # records, and two of the reference records are not connected to anything

        # in our first iteration of ref_to_ref, we should see the following pairs:
        # (ref0, ref2)
        # (ref0, ref3)
        # (ref1, ref2)
        # (ref1, ref3)

        # this is rough, so let's type it all out
        expected_pairs = set(
            [
                (ref_gbks[0].region, ref_gbks[2].region),
                (ref_gbks[0].region, ref_gbks[3].region),
                (ref_gbks[1].region, ref_gbks[2].region),
                (ref_gbks[1].region, ref_gbks[3].region),
            ]
        )

        actual_pairs = set(list(ref_to_ref_pair_generator.generate_pairs()))

        self.assertEqual(expected_pairs, actual_pairs)

    def test_ref_to_ref_pair_generator_second_iteration(self):
        """Tests whether the RefTorefPairGenerator correctly generates a set of
        pairs in the second iteration of a specific network
        """
        bs_data.DB.create_in_mem()

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]

        ref_to_ref_pair_generator = RefToRefRecordPairGenerator("mix", 1, "mix")
        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        # the state at the network at this point should be that all query to ref pairs
        # have been generated and their distances calculated. they are then added to
        # the network

        # the idea is that in the first iteration, this should compare all connected
        # ref nodes to all singleton ref nodes

        ref_to_ref_pair_generator.add_records(source_records)

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
                save_edge_to_db(
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
                save_edge_to_db(
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

        # so now we have a network where the query is connected to 2 of the reference
        # records, and two of the reference records are not connected to anything
        # let's do the first iteration
        list(ref_to_ref_pair_generator.generate_pairs())

        # we throw away the result because I want to test the second iteration, and
        # I want to enter the distance data manually

        # in our first iteration of ref_to_ref, we created the following pairs:
        # (ref0, ref2)
        # (ref0, ref3)
        # (ref1, ref2)
        # (ref1, ref3)

        # we will connect ref1 and ref3
        save_edge_to_db(
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
        save_edge_to_db(
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
        save_edge_to_db(
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
        save_edge_to_db(
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

        # so we compared query to all, and we compared ref0 to ref2 and ref3, and ref1
        # to ref2 and ref3. Each step has been an expansion of the previous step
        # now the next iteration should be to take the newly connected noe (ref3) and
        # compare it to all the other singletons (ref2)
        # so we only expect to see one pair:
        # (ref2, ref3)

        # now we can do the second iteration
        expected_pairs = set(
            [
                (ref_gbks[3].region, ref_gbks[2].region),
            ]
        )

        actual_pairs = set(list(ref_to_ref_pair_generator.generate_pairs()))

        self.assertEqual(expected_pairs, actual_pairs)

    def test_connected_component_pair_generator(self):
        """Tests whether the ConnectedComponenetPairGenerator correctly generates a set of
        pairs and memebers when given a connected component that only features a subset
        of the total members
        """

        bs_data.DB.create_in_mem()

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }
        weights = "mix"

        edge_param_id = bs_comparison.get_edge_param_id(run, weights)
        list([edge_param_id])

        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        # query -> test_path_0.gbk, rec_id 1
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]
        # ref[0] -> test_path_1.gbk, rec_id 2
        # ref[1] -> test_path_2.gbk, rec_id 3

        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        # making query <-> ref_1 edge with distance 0.0
        connected_component = [
            (
                query_gbk.region._db_id,
                ref_gbks[0].region._db_id,
                0.0,
                1.0,
                1.0,
                1.0,
                1,
            ),
            (
                query_gbk.region._db_id,
                ref_gbks[1].region._db_id,
                1.0,
                0.0,
                0.0,
                0.0,
                1,
            ),
            (
                ref_gbks[0].region._db_id,
                ref_gbks[1].region._db_id,
                0.0,
                1.0,
                1.0,
                1.0,
                1,
            ),
        ]

        expected_pairs = set(
            [
                (query_gbk.region, ref_gbks[0].region),
                (query_gbk.region, ref_gbks[1].region),
                (ref_gbks[0].region, ref_gbks[1].region),
            ]
        )
        # expected_record_ids = [1, 2, 3]

        cc_pair_generator = ConnectedComponentPairGenerator(connected_component, "mix")
        cc_pair_generator.add_records(source_records)

        # actual_record_ids = cc_pair_generator.record_ids = [1, 2, 3]
        actual_pairs = set(list(cc_pair_generator.generate_pairs()))

        self.assertEqual(expected_pairs, actual_pairs)

    def test_cull_singletons_cutoff(self):
        """Tests whether singletons are correctly culled"""

        bs_data.DB.create_in_mem()
        query_gbk = create_mock_gbk(0, bs_enums.SOURCE_TYPE.QUERY)
        # query -> test_path_0.gbk, rec_id 1
        query_gbk.save_all()
        ref_gbks = [
            create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE) for i in range(1, 5)
        ]
        # ref[0] -> test_path_1.gbk, rec_id 2
        # ref[1] -> test_path_2.gbk, rec_id 3

        source_records = [query_gbk.region]
        for ref_gbk in ref_gbks:
            source_records.append(ref_gbk.region)
            ref_gbk.save_all()

        new_bin = RecordPairGenerator("Test", weights="mix", edge_param_id=1)
        new_bin.add_records(source_records)

        # making query <-> ref_1 edge with distance 0.0

        save_edge_to_db(
            (
                query_gbk.region._db_id,
                ref_gbks[0].region._db_id,
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

        save_edge_to_db(
            (
                query_gbk.region._db_id,
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

        save_edge_to_db(
            (
                ref_gbks[0].region._db_id,
                ref_gbks[1].region._db_id,
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

        # edges above cutoff:
        # query <-> ref_1 | rec_id 1 <-> rec_id 2
        # ref_1 <-> ref_2 | rec_id 2 <-> rec_id 3

        new_bin.cull_singletons(0.5)

        expected_records = [source_records[0], source_records[1], source_records[2]]
        # expected_record_ids = set([1, 2, 3])

        actual_records = new_bin.source_records
        # actual_record_ids = new_bin.record_ids

        self.assertEqual(expected_records, actual_records)


class TestMixComparison(TestCase):
    def test_mix_iter(self):
        """Tests whether a new mix bin can be created for comparison"""
        gbk1 = GBK(Path("test"), "test1", source_type=bs_enums.SOURCE_TYPE.QUERY)
        gbk2 = GBK(Path("test"), "test2", source_type=bs_enums.SOURCE_TYPE.QUERY)
        gbk3 = GBK(Path("test"), "test3", source_type=bs_enums.SOURCE_TYPE.QUERY)

        bgc_a = BGCRecord(gbk1, 0, 0, 10, False, "")
        bgc_a.parent_gbk = gbk1
        bgc_a._db_id = 1

        bgc_b = BGCRecord(gbk2, 0, 0, 10, False, "")
        bgc_b.parent_gbk = gbk2
        bgc_b._db_id = 2

        bgc_c = BGCRecord(gbk3, 0, 0, 10, False, "")
        bgc_c.parent_gbk = gbk3
        bgc_c._db_id = 3

        bgc_list = [bgc_a, bgc_b, bgc_c]

        new_bin = generate_mix_bin(bgc_list, 1, bs_enums.RECORD_TYPE.REGION)

        # expected representation of the bin object
        expected_pair_count = 3

        actual_pair_count = len(list(new_bin.generate_pairs()))

        self.assertEqual(expected_pair_count, actual_pair_count)


def mock_region() -> Region:
    """generates a mock region for testing"""

    region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
    region_feature.qualifiers = {
        "region_number": ["1"],
        "candidate_cluster_numbers": ["1"],
        "product": ["T1PKS"],
    }

    region = Region.parse_as5(region_feature)

    candidate_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
    candidate_cluster_feature.qualifiers = {
        "candidate_cluster_number": ["1"],
        "kind": ["neighbouring"],
        "protoclusters": ["1"],
        "product": ["T1PKS"],
    }

    candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

    protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature.qualifiers = {
        "protocluster_number": ["1"],
        "category": ["PKS"],
        "product": ["T1PKS"],
    }

    protocluster = ProtoCluster.parse(protocluster_feature)

    protocore_feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature.qualifiers = {
        "protocluster_number": ["1"],
        "product": ["T1PKS"],
    }

    protocore = ProtoCore.parse(protocore_feature)

    protocluster.add_proto_core(protocore)
    candidate_cluster.add_proto_cluster(protocluster)
    region.add_cand_cluster(candidate_cluster)

    return region


class TestBinGenerators(TestCase):
    """Test class for the bin generators and associated functions"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_get_region_category(self):
        """Tests whether a category is correclty parsed from a region of version as6 or higher"""

        region = mock_region()
        cc = region.cand_clusters[1]
        pc = cc.proto_clusters[1]

        expected_category = "PKS"
        category = get_record_category(pc)

        self.assertEqual(expected_category, category)

    def test_get_legacy_weight_from_category(self):
        """Tests wether the correct legacy weight category is created from a region category"""

        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CLASS,
            "record_type": bs_enums.RECORD_TYPE.REGION,
            "hybrids_off": False,
        }

        region = mock_region()
        cc = region.cand_clusters[1]
        pc = cc.proto_clusters[1]

        expected_category = "T1PKS"
        category = get_legacy_weights_from_category(pc, "T1PKS", run)

        self.assertEqual(expected_category, category)

    def test_as_class_bin_generator(self):
        """Tests whether an antismash bin is correclty generated given a weight label"""

        "Bin 'PKS': 1 pairs from 2 BGC records"
        bs_data.DB.create_in_mem()

        gbk_1 = GBK(Path("test"), "test", source_type=bs_enums.SOURCE_TYPE.QUERY)
        gbk_1.region = mock_region()

        region_1 = gbk_1.region
        cand_cluster_1 = region_1.cand_clusters[1]
        protocluster_1 = cand_cluster_1.proto_clusters[1]
        protocore_1 = protocluster_1.proto_core[1]

        gbk_2 = GBK(Path("test"), "test", source_type=bs_enums.SOURCE_TYPE.QUERY)
        gbk_2.region = mock_region()

        region_2 = gbk_2.region
        cand_cluster_2 = region_2.cand_clusters[1]
        protocluster_2 = cand_cluster_2.proto_clusters[1]
        protocore_2 = protocluster_2.proto_core[1]

        # gbks = [gbk_1.region, gbk_2.region]
        gbks = [protocore_1, protocore_2]

        expected_dict = OrderedDict()
        expected_dict = {
            "classmix": "T1PKSmix",
            "categorymix": "PKSmix",
            "classlegacy_weights": "T1PKST1PKS",
            "categorylegacy_weights": "PKST1PKS",
        }

        # run_class_mix = {
        #     "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
        #     "legacy_weights": False,
        #     "classify": bs_enums.CLASSIFY_MODE.CLASS,
        # }

        run_category_weights = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
            "hybrids_off": False,
            "record_type": bs_enums.RECORD_TYPE.REGION,
        }

        bin = next(as_class_bin_generator(gbks, run_category_weights))
        expected_combo = expected_dict["categorylegacy_weights"]
        seen_combo = bin.label + bin.weights

        self.assertEqual(expected_combo, seen_combo)

    def test_get_edge_params_id_insert(self):
        """Tests insertion of edge params into the database when not there"""

        bs_data.DB.create_in_mem()
        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }
        weights = "mix"

        expected_id = 1
        actual_id = bs_comparison.get_edge_param_id(run, weights)

        self.assertEqual(expected_id, actual_id)

    def test_get_edge_params_id_fetch(self):
        """Tests getting edge params from the database when already there"""

        bs_data.DB.create_in_mem()
        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }

        first_id = bs_comparison.get_edge_param_id(run, "mix")
        list([first_id])
        second_id = bs_comparison.get_edge_param_id(run, "mix")

        expected_id = 1

        self.assertEqual(expected_id, second_id)

    def test_get_edge_weight(self):
        """Tests whether an edge weight gets correctly fetched from the database"""

        bs_data.DB.create_in_mem()
        run = {
            "alignment_mode": bs_enums.ALIGNMENT_MODE.AUTO,
            "legacy_weights": True,
            "classify": bs_enums.CLASSIFY_MODE.CATEGORY,
        }

        first_id = bs_comparison.get_edge_param_id(run, "mix")
        second_id = bs_comparison.get_edge_param_id(run, "other")
        list([first_id, second_id])

        weights = bs_comparison.get_edge_weight(second_id)
        expected_weights = "other"

        self.assertEqual(expected_weights, weights)

    def test_get_query_records(self):
        self.skipTest("Not implemented")
