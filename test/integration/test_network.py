"""module containing integration tests for the distance generatio/comparison workflows"""

# from python
from unittest import TestCase
from pathlib import Path
from itertools import combinations

from sqlalchemy import select

# from other modules
import big_scape.data as bs_data
import big_scape.comparison as bs_comparison
import big_scape.genbank as bs_gbk
import big_scape.network.network as bs_network
from big_scape.genbank import (
    GBK,
    CDS,
    Region,
)
import big_scape.enums as bs_enums
import big_scape.file_input as bs_files


def create_mock_gbk(i, source_type: bs_enums.SOURCE_TYPE) -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, i, 0, 100, False, "test")
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
    score: float = 0.0,
) -> list[
    tuple[int, int, float, float, float, float, int, bs_comparison.ComparableRegion]
]:
    edges = []
    for gbk_a, gbk_b in combinations(edge_gbks, 2):
        if gbk_a.region is None or gbk_b.region is None:
            continue
        if gbk_a.region._db_id is None or gbk_b.region._db_id is None:
            continue

        edges.append(
            edge_creation_function(gbk_a.region._db_id, gbk_b.region._db_id, score)
        )

    return edges


def create_mock_edge(a_id, b_id, distance=0.0):
    return (
        a_id,
        b_id,
        distance,
        1 - distance,
        1 - distance,
        1 - distance,
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
    """Test class for the network module"""

    def clean_db(self):
        if bs_data.DB.opened():
            bs_data.DB.close_db()
        bs_data.DB.metadata = None

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)
        self.addCleanup(self.clean_db)

    def test_get_connected_components_two_cutoffs(self):
        bs_data.DB.create_in_mem()

        gbks_a = []
        for i in range(3):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_a.append(gbk)
            gbk.save_all()

        gbks_b = []
        for i in range(3):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_b.append(gbk)
            gbk.save_all()

        # create a bunch of edges, generating a cc with edges of distance 0
        # and another with distance 0.6
        # gen_mock_edge_list creates edges for all combinations of gbks
        edges_a = gen_mock_edge_list(gbks_a, create_mock_edge)
        edges_b = gen_mock_edge_list(gbks_b, create_mock_edge, 0.6)

        edges = edges_a + edges_b

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        bs_network.get_connected_components(0.5, 1, gbks_a + gbks_b)

        ccs = list(bs_network.get_connected_components(0.5, 1, gbks_a + gbks_b))

        self.assertEqual(len(ccs), 1)

        ccs = list(bs_network.get_connected_components(1, 1, gbks_a + gbks_b))

        self.assertEqual(len(ccs), 2)

    def test_get_connected_components_all_cutoffs(self):
        bs_data.DB.create_in_mem()

        gbks_a = []
        for i in range(3):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_a.append(gbk)
            gbk.save_all()

        gbks_b = []
        for i in range(3, 6):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_b.append(gbk)
            gbk.save_all()

        gbks_c = []
        for i in range(6, 9):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_c.append(gbk)
            gbk.save_all()

        all_gbks = gbks_a + gbks_b + gbks_c

        # create a bunch of edges, generating a cc with edges of distance 0
        # and another with distance 0.6
        # gen_mock_edge_list creates edges for all combinations of gbks

        edges_a = gen_mock_edge_list(gbks_a, create_mock_edge)
        edges_b = gen_mock_edge_list(gbks_b, create_mock_edge, 0.6)
        edges_c = gen_mock_edge_list(gbks_c, create_mock_edge, 0.2)

        edges_d = gen_mock_edge_list(
            [gbks_a[0], gbks_b[0], gbks_c[0]], create_mock_edge, 0.8
        )

        edges = edges_a + edges_b + edges_c + edges_d

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # distance cutoff 0.0
        ccs = list(bs_network.get_connected_components(0.0, 1, all_gbks))

        self.assertEqual(len(ccs), 0)

        # distance cutoff 0.1
        ccs = list(bs_network.get_connected_components(0.1, 1, all_gbks))

        self.assertEqual(len(ccs), 1)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.1)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, len(gbks_a))
        self.assertEqual(total_edges_in_ccs, len(edges_a))

        # distance cutoff 0.2
        ccs = list(bs_network.get_connected_components(0.2, 1, all_gbks))

        self.assertEqual(len(ccs), 1)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.2)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, len(gbks_a))
        self.assertEqual(total_edges_in_ccs, len(edges_a))

        # distance cutoff 0.3
        ccs = list(bs_network.get_connected_components(0.3, 1, all_gbks))

        self.assertEqual(len(ccs), 2)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.3)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, len(gbks_a + gbks_b))
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b))

        # distance cutoff 0.4
        ccs = list(bs_network.get_connected_components(0.4, 1, all_gbks))

        self.assertEqual(len(ccs), 2)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.4)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, len(gbks_a + gbks_b))
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b))

        # distance cutoff 0.5
        ccs = list(bs_network.get_connected_components(0.5, 1, all_gbks))

        self.assertEqual(len(ccs), 2)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.5)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, 6)
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b))

        # distance cutoff 0.6
        ccs = list(bs_network.get_connected_components(0.6, 1, all_gbks))

        self.assertEqual(len(ccs), 2)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.6)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, 6)
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b))

        # distance cutoff 0.7
        ccs = list(bs_network.get_connected_components(0.7, 1, all_gbks))

        self.assertEqual(len(ccs), 3)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.7)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, 9)
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b + edges_c))

        # distance cutoff 0.8
        ccs = list(bs_network.get_connected_components(0.8, 1, all_gbks))

        self.assertEqual(len(ccs), 3)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.8)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, 9)
        self.assertEqual(total_edges_in_ccs, len(edges_a + edges_b + edges_c))

        # distance cutoff 0.9
        ccs = list(bs_network.get_connected_components(0.9, 1, all_gbks))

        self.assertEqual(len(ccs), 1)

        total_nodes_in_ccs = 0
        total_edges_in_ccs = 0

        for cc in ccs:
            nodes_in_cc = len(bs_network.get_nodes_from_cc(cc, all_gbks))
            total_nodes_in_ccs += nodes_in_cc

            for edge in cc:
                self.assertLess(edge[2], 0.9)
                total_edges_in_ccs += 1

        self.assertEqual(total_nodes_in_ccs, 9)
        self.assertEqual(total_edges_in_ccs, len(edges))

        # check state of the cctable at this point, which contains one row
        # per node per connected component per cutoff

        cc_table = bs_data.DB.metadata.tables["connected_component"]

        select_statement = select(cc_table.c.id, cc_table.c.cutoff)

        cc_ids_and_cutoffs = bs_data.DB.execute(select_statement).fetchall()

        len(cc_ids_and_cutoffs)

        expected_len = (
            9 * 0
            + 3 * 1  # distance 0.0, 9 nodes in 0 ccs
            + 3 * 1  # distance 0.1, 3 nodes in 1 cc
            + 3 * 2  # distance 0.2, 3 nodes in 1 cc
            + 3 * 2  # distance 0.3, 3 nodes in 2 ccs
            + 3 * 2  # distance 0.4, 3 nodes in 2 ccs
            + 3 * 2  # distance 0.5, 3 nodes in 2 ccs
            + 3 * 3  # distance 0.6, 3 nodes in 2 ccs
            + 3 * 3  # distance 0.7, 3 nodes in 3 ccs
            + 9 * 1  # distance 0.8, 3 nodes in 3 ccs  # distance 0.9, 9 nodes in 1 cc
        )

        self.assertEqual(len(cc_ids_and_cutoffs), expected_len)

    def test_get_connected_components_no_ref_to_ref_ccs(self):
        """Tests whether ref only ccs are correclty excluded from the analysis"""
        bs_data.DB.create_in_mem()
        run = {"record_type": bs_enums.RECORD_TYPE.REGION}

        gbks_a = []
        for i in range(3):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE)
            gbks_a.append(gbk)
            gbk.save_all()

        gbks_b = []
        for i in range(3, 6):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.QUERY)
            gbks_b.append(gbk)
            gbk.save_all()

        gbks_c = []
        for i in range(6, 9):
            gbk = create_mock_gbk(i, bs_enums.SOURCE_TYPE.REFERENCE)
            gbks_c.append(gbk)
            gbk.save_all()

        all_gbks = gbks_a + gbks_b + gbks_c

        include_records = []

        record = bs_files.get_all_bgc_records(run, all_gbks)
        include_records.extend(record)

        # create a bunch of edges, generating a cc with edges of distance 0
        # and another with distance 0.6
        # gen_mock_edge_list creates edges for all combinations of gbks

        edges_a = gen_mock_edge_list(gbks_a, create_mock_edge)
        edges_b = gen_mock_edge_list(gbks_b, create_mock_edge)
        edges_c = gen_mock_edge_list(gbks_c, create_mock_edge, 0.6)

        edges_d = gen_mock_edge_list([gbks_b[0], gbks_c[0]], create_mock_edge, 0.2)

        edges = edges_a + edges_b + edges_c + edges_d

        # save the edges
        for edge in edges:
            bs_comparison.save_edge_to_db(edge)

        # distance cutoff 0.8, here we should only get one connected component,
        # featuring all nodes from gbks_b and gbks_c
        ccs = bs_network.get_connected_components(0.8, 1, all_gbks)

        for cc in ccs:
            is_ref_only = bs_network.reference_only_connected_component(
                cc, include_records
            )
            if is_ref_only:
                bs_network.remove_connected_component(cc, 0.8, 1)

        cc_table = bs_data.DB.metadata.tables["connected_component"]

        select_statement = select(cc_table.c.id).distinct()
        len_cc = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(len_cc, 1)

        select_statement = select(cc_table.c.record_id).where(cc_table.c.cutoff == 0.8)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, len(gbks_b + gbks_c))

        # distance cutoff 0.5, here we should only get one connected component,
        # featuring nodes from gbks_b and one node from gbks_c
        ccs = bs_network.get_connected_components(0.5, 1, all_gbks)

        for cc in ccs:
            is_ref_only = bs_network.reference_only_connected_component(
                cc, include_records
            )
            if is_ref_only:
                bs_network.remove_connected_component(cc, 0.5, 1)

        cc_table = bs_data.DB.metadata.tables["connected_component"]

        select_statement = select(cc_table.c.id).distinct()
        len_cc = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(len_cc, 1)

        select_statement = select(cc_table.c.record_id).where(cc_table.c.cutoff == 0.5)
        rows = len(bs_data.DB.execute(select_statement).fetchall())

        self.assertEqual(rows, len(gbks_b + [gbks_c[0]]))
