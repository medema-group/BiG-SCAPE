"""Contains tests that check if the partial data state functions correctly determine and
return the state of partial analyses so that they can be continued
"""

# from python
from pathlib import Path
from unittest import TestCase
from itertools import combinations

# from other modules
from big_scape.data import (
    DB,
    find_minimum_task,
    get_input_data_state,
    get_hmm_data_state,
    get_comparison_data_state,
    get_cds_to_scan,
    get_hsp_to_align,
)
from big_scape.genbank import GBK, CDS, Region
from big_scape.hmm import HSP, HSPAlignment, HMMer

import big_scape.comparison as bs_comparison

import big_scape.enums as bs_enums


def create_mock_gbk(i) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), bs_enums.SOURCE_TYPE.QUERY)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
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


def add_mock_hsp_cds(cds: CDS) -> None:
    hsp = HSP(cds, "PF01234", 1.0, 0, 100)
    cds.hsps.append(hsp)


def add_mock_hsp_alignment_hsp(hsp: HSP) -> None:
    hsp_alignment = HSPAlignment(hsp, "AAAAAAAAAA")
    hsp.alignment = hsp_alignment


def gen_mock_edge_list(
    edge_gbks: list[GBK],
) -> list[
    tuple[
        int,
        int,
        float,
        float,
        float,
        float,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        int,
        bool,
    ]
]:
    edges = []
    for gbk_a, gbk_b in combinations(edge_gbks, 2):
        if gbk_a.region is None or gbk_b.region is None:
            continue
        if gbk_a.region._db_id is None or gbk_b.region._db_id is None:
            continue

        edges.append(
            (
                gbk_a.region._db_id,
                gbk_b.region._db_id,
                0.0,
                0.0,
                0.0,
                0.0,
                1,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                False,
            )
        )

    return edges


class TestPartial(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_min_task_data(self):
        DB.create_in_mem()
        expected_min_task = bs_enums.TASK.SAVE_GBKS

        gbks = [create_mock_gbk(1)]

        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_scan(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(1)]
        # add gbk to db
        gbks[0].save_all()

        expected_min_task = bs_enums.TASK.HMM_SCAN
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_hmm_align(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(1)]
        # add hsp to gbk
        add_mock_hsp_cds(gbks[0].genes[0])

        # add gbk to db
        gbks[0].save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])

        expected_min_task = bs_enums.TASK.HMM_ALIGN
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)

    def test_min_task_comparison(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(1)]
        # add hsp to gbk
        add_mock_hsp_cds(gbks[0].genes[0])
        # add hspalignment to hsp
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])

        # add gbk to db
        gbks[0].save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])
        # sigh
        gbks[0].genes[0].hsps[0].alignment.save()

        expected_min_task = bs_enums.TASK.COMPARISON
        actual_min_task = find_minimum_task(gbks)

        self.assertEqual(expected_min_task, actual_min_task)


class TestPartialInputs(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_data(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(1)]

        expected_state = bs_enums.INPUT_TASK.NO_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_partial_data(self):
        DB.create_in_mem()

        # gbk 1
        gbk_1_in_db = create_mock_gbk(1)
        gbk_1_in_db.save_all()
        # gbk 2
        gbk_2_in_db = create_mock_gbk(2)
        gbk_2_in_db.save_all()

        gbks = [gbk_1_in_db]

        expected_state = bs_enums.INPUT_TASK.PARTIAL_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_data(self):
        DB.create_in_mem()

        gbk_in_db = create_mock_gbk(1)
        gbk_in_db.save_all()

        gbk_not_in_db = create_mock_gbk(2)
        gbks = [gbk_in_db, gbk_not_in_db]

        expected_state = bs_enums.INPUT_TASK.MIXED_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_no_new_data(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
            create_mock_gbk(3),
        ]
        for gbk in gbks:
            gbk.save_all()

        expected_state = bs_enums.INPUT_TASK.SAME_DATA
        actual_state = get_input_data_state(gbks)

        self.assertEqual(expected_state, actual_state)


class TestPartialHMM(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_scans_done(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(1)]

        expected_state = bs_enums.HMM_TASK.NO_DATA
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_scans_to_do(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
        ]
        # add hsp to one gbk
        add_mock_hsp_cds(gbks[0].genes[0])

        # add gbks to db
        for gbk in gbks:
            gbk.save_all()
        # add hsp to db
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])

        expected_state = bs_enums.HMM_TASK.NEED_SCAN
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_align_to_do(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
        ]
        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        # add hspalignment to only one hsp
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])

        # add gbks and hsps to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
        # add hspalignment to db
        gbks[0].genes[0].hsps[0].alignment.save()

        expected_state = bs_enums.HMM_TASK.NEED_ALIGN
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_no_new_align(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
        ]
        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        # add hspalignment to both hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        expected_state = bs_enums.HMM_TASK.ALL_ALIGNED
        actual_state = get_hmm_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_get_cds_to_scan(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
        ]
        # add hsp to first gbk
        add_mock_hsp_cds(gbks[0].genes[0])

        # add gbks and cds to db
        for gbk in gbks:
            gbk.save_all()

        # save first hsp and only set the first cds as scanned
        gbks[0].genes[0].hsps[0].save()
        HMMer.set_hmm_scanned(gbks[0].genes[0])

        expected_cds_to_scan = [gbks[1].genes[0]]
        actual_cds_to_scan = get_cds_to_scan(gbks)

        self.assertEqual(expected_cds_to_scan, actual_cds_to_scan)

    def test_get_hsp_to_align(self):
        DB.create_in_mem()

        gbks = [
            create_mock_gbk(1),
            create_mock_gbk(2),
        ]
        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])

        # add gbks and cds to db. set cds as scnaned
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbks[0].genes[0])

        # add alignment only to first hsp
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])

        expected_hsp_to_align = [gbks[1].genes[0].hsps[0]]
        actual_hsp_to_align = get_hsp_to_align(gbks)

        self.assertEqual(expected_hsp_to_align, actual_hsp_to_align)


class TestPartialComparison(TestCase):
    def clean_db(self):
        if DB.opened():
            DB.close_db()
        DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_no_comparisons_done(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(i) for i in range(3)]

        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        add_mock_hsp_cds(gbks[2].genes[0])
        # add hspalignment to both hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[2].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        # no alignments done. expect no data

        expected_state = bs_enums.COMPARISON_TASK.NO_DATA
        actual_state = get_comparison_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_new_comparisons_to_do(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(i) for i in range(3)]

        # add hsp to both gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        add_mock_hsp_cds(gbks[2].genes[0])
        # add hspalignment to both hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[2].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        # only one distance done. (1-2, missing 1-3 and 2-3)
        edges = gen_mock_edge_list(gbks[0:2])
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        expected_state = bs_enums.COMPARISON_TASK.NEW_DATA
        actual_state = get_comparison_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_no_new_comparisons(self):
        self.skipTest("Broken")
        DB.create_in_mem()

        gbks = [create_mock_gbk(i) for i in range(3)]

        # add hsp to gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        add_mock_hsp_cds(gbks[2].genes[0])
        # add hspalignment to hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[2].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        # only one distance done. (1-2, missing 1-3 and 2-3)
        edges = gen_mock_edge_list(gbks[0:3])
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        expected_state = bs_enums.COMPARISON_TASK.ALL_DONE
        actual_state = get_comparison_data_state(gbks)

        self.assertEqual(expected_state, actual_state)

    def test_partial_pair_generator(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(i) for i in range(3)]

        # add hsp to gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        add_mock_hsp_cds(gbks[2].genes[0])
        # add hspalignment to hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[2].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        # only one distance done. (1-2, missing 1-3 and 2-3)
        edges = gen_mock_edge_list(gbks[0:2])
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # all-vs-all bin
        mix_bin = bs_comparison.RecordPairGenerator("mix", 1)
        mix_bin.add_records([gbk.region for gbk in gbks])

        expected_missing_pairs = [
            bs_comparison.RecordPair(gbks[0].region, gbks[2].region),
            bs_comparison.RecordPair(gbks[1].region, gbks[2].region),
        ]

        pair_generator = bs_comparison.RecordPairGenerator("mix", 1)
        pair_generator.add_records([gbk.region for gbk in gbks])

        missing_edge_generator = bs_comparison.MissingRecordPairGenerator(
            pair_generator
        )

        actual_missing_pairs = list(missing_edge_generator.generate_pairs())

        self.assertListEqual(expected_missing_pairs, actual_missing_pairs)

    def test_get_missing_distance_count(self):
        DB.create_in_mem()

        gbks = [create_mock_gbk(i) for i in range(3)]

        # add hsp to gbks
        add_mock_hsp_cds(gbks[0].genes[0])
        add_mock_hsp_cds(gbks[1].genes[0])
        add_mock_hsp_cds(gbks[2].genes[0])
        # add hspalignment to hsps
        add_mock_hsp_alignment_hsp(gbks[0].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[1].genes[0].hsps[0])
        add_mock_hsp_alignment_hsp(gbks[2].genes[0].hsps[0])

        # add gbks, hsps and alignments to db
        for gbk in gbks:
            gbk.save_all()
            gbk.genes[0].hsps[0].save()
            HMMer.set_hmm_scanned(gbk.genes[0])
            gbk.genes[0].hsps[0].alignment.save()

        # only one distance done. (1-2, missing 1-3 and 2-3)
        edges = gen_mock_edge_list(gbks[0:2])
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # all-vs-all bin
        mix_bin = bs_comparison.RecordPairGenerator("mix", 1)
        mix_bin.add_records([gbk.region for gbk in gbks])

        expected_missing_count = 2

        pair_generator = bs_comparison.RecordPairGenerator("mix", 1)
        pair_generator.add_records([gbk.region for gbk in gbks])

        missing_edge_generator = bs_comparison.MissingRecordPairGenerator(
            pair_generator
        )

        actual_missing_count = missing_edge_generator.num_pairs()

        self.assertEqual(expected_missing_count, actual_missing_count)
