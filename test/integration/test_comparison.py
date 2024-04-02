"""module containing integration tests for the distance generatio/comparison workflows"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
import big_scape.data as bs_data
import big_scape.comparison as bs_comparison
from big_scape.genbank import GBK, CDS, Region
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

    def test_calculate_distances_legacy_classify(self):
        pass

    def test_calculate_distances_classify(self):
        pass

    def test_calculate_distances_query(self):
        pass
