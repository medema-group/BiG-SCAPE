"""Contains tests for extension of comparable regions with record pairs"""

# from python
import unittest
import random
import string

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.hmm as bs_hmmer
import big_scape.comparison as bs_comp
import big_scape.enums as bs_enums


def generate_random_hsp(cds: bs_genbank.CDS):
    """Generate a random domain for a cds"""
    # generate a random string as accession
    accession = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    hsp = bs_hmmer.HSP(cds, accession, 100.0, cds.nt_start, cds.nt_stop)
    return hsp


def generate_common_hsp(cds_a: bs_genbank.CDS, cds_b: bs_genbank.CDS, offset=0):
    """Generate a random domain common between a pair of cds"""
    # generate a random string as accession
    accession = "".join(random.choices(string.ascii_uppercase + string.digits, k=10))
    hsp_a = bs_hmmer.HSP(
        cds_a, accession, 100.0, cds_a.nt_start + offset, cds_a.nt_stop + offset
    )
    hsp_b = bs_hmmer.HSP(
        cds_b, accession, 100.0, cds_b.nt_start + offset, cds_b.nt_stop + offset
    )

    return hsp_a, hsp_b


def generate_mock_cds_lists(
    len_a: int, len_b: int, common_a: list[int], common_b: list[int], reverse: bool
):
    """Generate a mock cds list"""

    if len(common_a) != len(common_b):
        raise ValueError("common_a and common_b must be of the same length")

    cds_a: list[bs_genbank.CDS] = []
    cds_b: list[bs_genbank.CDS] = []
    for i in range(max(len_a, len_b)):
        if i < len_a:
            a = bs_genbank.CDS(i * 100, (i + 1) * 100)
            cds_a.append(a)
        if i < len_b:
            b = bs_genbank.CDS(i * 100, (i + 1) * 100)
            cds_b.append(b)

        if i not in common_a:
            cds_a[i].hsps = [generate_random_hsp(cds_a[i])]

        if i not in common_b:
            cds_b[i].hsps = [generate_random_hsp(cds_b[i])]

    for common_idx in range(len(common_a)):
        common_a_idx = common_a[common_idx]
        common_b_idx = common_b[common_idx]

        common_a_hsp, common_b_hsp = generate_common_hsp(
            cds_a[common_a_idx], cds_b[common_b_idx], common_idx * 10
        )
        cds_a[common_a_idx].hsps.append(common_a_hsp)
        cds_b[common_b_idx].hsps.append(common_b_hsp)

    if reverse:
        cds_b = cds_b[::-1]
        for cds in cds_b:
            cds.hsps = cds.hsps[::-1]

    return cds_a, cds_b


def generate_mock_region(cds_list: list[bs_genbank.CDS]) -> bs_genbank.Region:
    """Generate a mock region from a cds list"""
    gbk = bs_genbank.GBK(None, "", bs_enums.SOURCE_TYPE.QUERY)
    gbk.genes = cds_list

    region = bs_genbank.Region(gbk, 1, 0, len(cds_list) * 100, False, "")

    return region


def generate_mock_protocluster(cds_list: list[bs_genbank.CDS], protocore_idx: int):
    """Generate a mock protocluster from a cds list"""
    gbk = bs_genbank.GBK(None, "", bs_enums.SOURCE_TYPE.QUERY)
    gbk.genes = cds_list
    protocluster = bs_genbank.ProtoCluster(
        gbk, 1, 0, len(cds_list) * 100, False, "", {}, ""
    )

    protocore_start = cds_list[protocore_idx].nt_start
    protocore_stop = cds_list[protocore_idx].nt_stop

    protocore = bs_genbank.ProtoCore(gbk, 1, protocore_start, protocore_stop, False, "")

    protocluster.proto_core[1] = None
    protocluster.add_proto_core(protocore)

    return protocluster


class TestExtendUtilities(unittest.TestCase):
    """Tests for extension utilities"""

    def test_reset_region(self):
        """Test for reset method on a region"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_region(cds_a)
        record_b = generate_mock_region(cds_b)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        # reset the expansion
        bs_comp.extend.reset(comparable_region)

        expected_comparable_region = bs_comp.ComparableRegion(
            record_pair, 0, len(cds_a), 0, len(cds_b), False
        )

        self.assertEqual(expected_comparable_region, comparable_region)

    def test_reset_protocluster(self):
        """Test for reset method on a protocluster"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        # reset the expansion
        bs_comp.extend.reset(comparable_region)

        expected_comparable_region = bs_comp.ComparableRegion(
            record_pair, 0, len(cds_a), 0, len(cds_b), False
        )

        self.assertEqual(expected_comparable_region, comparable_region)

    def test_check_too_short(self):
        """Test for check method when a region is too short and contains no
        biosynthetic genes
        """

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        expected_result = False
        actual_result = bs_comp.extend.check(comparable_region, 3, True)

        self.assertEqual(expected_result, actual_result)

    def test_check_too_short_biosynth(self):
        """Test for check method when the region contains a biosynthetic cds"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        cds_a[5].gene_kind = "biosynthetic"

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        expected_result = True
        actual_result = bs_comp.extend.check(comparable_region, 3, True)

        self.assertEqual(expected_result, actual_result)

    def test_check_pass(self):
        """Test for check method when the region is long enough without a biosynthetic
        cds
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 9, 5, 9, False)

        expected_result = True
        actual_result = bs_comp.extend.check(comparable_region, 3, True)

        self.assertEqual(expected_result, actual_result)


class TestScoreExtend(unittest.TestCase):
    """Tests for score extension"""

    def test_get_target_indexes(self):
        """Tests for getting target domain (cds) indexes"""

        # create mock cds list
        target = []
        for i in range(5):
            cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
            target.append(cds)

        # add mock domains
        target[0].hsps = [bs_hmmer.HSP(target[0], "A", 100.0, 0, 100)]
        target[1].hsps = [bs_hmmer.HSP(target[1], "B", 100.0, 0, 100)]
        target[2].hsps = [bs_hmmer.HSP(target[2], "C", 100.0, 0, 100)]
        target[3].hsps = [bs_hmmer.HSP(target[3], "B", 100.0, 0, 100)]
        # add a duplicate to the same cds
        target[3].hsps.append(bs_hmmer.HSP(target[3], "B", 100.0, 0, 100))
        target[4].hsps = [bs_hmmer.HSP(target[4], "A", 100.0, 0, 100)]
        # add a duplicate to a different cds
        target[4].hsps.append(bs_hmmer.HSP(target[4], "C", 100.0, 0, 100))

        expected_target_index = {
            "A": [(0, 0), (4, 5)],
            "B": [(1, 1), (3, 3), (3, 4)],
            "C": [(2, 2), (4, 6)],
        }

        actual_target_index = bs_comp.extend.get_target_indexes(target)

        self.assertEqual(expected_target_index, actual_target_index)

    def test_extend_all_match(self):
        """Tests for extension with all matches"""
        # all cds have a common domains
        query, target = generate_mock_cds_lists(
            5, 5, [0, 1, 2, 3, 4], [0, 1, 2, 3, 4], False
        )

        target_index = bs_comp.extend.get_target_indexes(target)

        expected_extends = (5, 5, 25)

        actual_extends = bs_comp.extend.score_extend(
            query, 0, target_index, 5, -5, -3, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_large_gap(self):
        """Tests for extension with a large gap"""
        # only first and last have common domains
        # so gap of 3. 3*-2 = -6 with one match = -1. so no extension
        query = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(400, 500),
        ]
        target = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(100, 200),
            bs_genbank.CDS(200, 300),
            bs_genbank.CDS(300, 400),
            bs_genbank.CDS(400, 500),
        ]

        # common domains
        query[0].hsps = [bs_hmmer.HSP(query[0], "A", 100.0, 0, 100)]
        target[0].hsps = [bs_hmmer.HSP(target[0], "A", 100.0, 0, 100)]
        query[1].hsps = [bs_hmmer.HSP(query[1], "B", 100.0, 0, 100)]
        target[4].hsps = [bs_hmmer.HSP(target[4], "B", 100.0, 0, 100)]

        # different domains. only for target
        target[1].hsps = [bs_hmmer.HSP(target[1], "C", 100.0, 0, 100)]
        target[2].hsps = [bs_hmmer.HSP(target[2], "D", 100.0, 0, 100)]
        target[3].hsps = [bs_hmmer.HSP(target[3], "E", 100.0, 0, 100)]

        target_index = bs_comp.extend.get_target_indexes(target)

        expected_extends = (1, 1, 5)

        actual_extends = bs_comp.extend.score_extend(
            query, 0, target_index, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_small_gap(self):
        """Tests for extension with a small gap"""
        # two gaps, 3 matches. 3 * 5 + 2 * -2 = 11. so extend the full range
        query = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(100, 200),
            bs_genbank.CDS(400, 500),
        ]
        target = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(100, 200),
            bs_genbank.CDS(200, 300),
            bs_genbank.CDS(300, 400),
            bs_genbank.CDS(400, 500),
        ]

        # common domains
        query[0].hsps = [bs_hmmer.HSP(query[0], "A", 100.0, 0, 100)]
        target[0].hsps = [bs_hmmer.HSP(target[0], "A", 100.0, 0, 100)]
        query[1].hsps = [bs_hmmer.HSP(query[1], "B", 100.0, 0, 100)]
        target[1].hsps = [bs_hmmer.HSP(target[1], "B", 100.0, 0, 100)]
        query[2].hsps = [bs_hmmer.HSP(query[1], "C", 100.0, 0, 100)]
        target[4].hsps = [bs_hmmer.HSP(target[4], "C", 100.0, 0, 100)]

        # different domains. only for target
        target[2].hsps = [bs_hmmer.HSP(target[2], "D", 100.0, 0, 100)]
        target[3].hsps = [bs_hmmer.HSP(target[3], "E", 100.0, 0, 100)]

        target_index = bs_comp.extend.get_target_indexes(target)

        expected_extends = (5, 3, 11)

        actual_extends = bs_comp.extend.score_extend(
            query, 0, target_index, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_fill_gap(self):
        """Tests the concept of filling a gap during extend"""

        # create mock cds lists
        query = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(100, 200),
            bs_genbank.CDS(200, 300),
        ]
        target = [
            bs_genbank.CDS(0, 100),
            bs_genbank.CDS(100, 200),
            bs_genbank.CDS(200, 300),
        ]

        # the idea is that we have a situation like this:
        # query:  A B C
        # target: A C B
        # when we see something like this, we will have a gap of 1 when we reach B in
        # the query:
        # query:  A - B C
        # target: A C B
        # but we want to "fix" it when we see that close by there is another same domain
        # that would fill the gap

        query[0].hsps = [bs_hmmer.HSP(query[0], "A", 100.0, 0, 100)]
        target[0].hsps = [bs_hmmer.HSP(target[0], "A", 100.0, 0, 100)]

        query[1].hsps = [bs_hmmer.HSP(query[1], "B", 100.0, 0, 100)]
        target[1].hsps = [bs_hmmer.HSP(target[1], "C", 100.0, 0, 100)]

        query[2].hsps = [bs_hmmer.HSP(query[2], "C", 100.0, 0, 100)]
        target[2].hsps = [bs_hmmer.HSP(target[2], "B", 100.0, 0, 100)]

        # so we expect an extension of 3 on both sides, and a score of 15
        expected_extends = (3, 3, 15)

        target_index = bs_comp.extend.get_target_indexes(target)

        actual_extends = bs_comp.extend.score_extend(
            query, 0, target_index, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)
