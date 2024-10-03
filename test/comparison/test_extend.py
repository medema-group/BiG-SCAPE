"""Contains tests for extension of comparable regions with record pairs"""

# from python
import unittest
import random
import string
import big_scape.comparison.record_pair

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.hmm as bs_hmmer
import big_scape.comparison as bs_comp
import big_scape.enums as bs_enums


def generate_random_hsp(cds: bs_genbank.CDS):
    """Generate a random domain for a cds"""
    # generate a random string as accession
    accession = "".join(random.choices(string.ascii_uppercase + string.digits, k=7))
    hsp = bs_hmmer.HSP(cds, accession, 100.0, cds.nt_start, cds.nt_stop)
    return hsp


def generate_common_hsp(cds_a: bs_genbank.CDS, cds_b: bs_genbank.CDS, offset=0):
    """Generate a random domain common between a pair of cds"""
    # generate a random string as accession
    accession = "".join(random.choices(string.ascii_uppercase + string.digits, k=7))
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
            a.strand = 1
            a.orf_num = i
            cds_a.append(a)
            if i not in common_a:
                cds_a[i].hsps = [generate_random_hsp(cds_a[i])]

        if i < len_b:
            b = bs_genbank.CDS(i * 100, (i + 1) * 100)
            b.strand = 1
            b.orf_num = i
            cds_b.append(b)
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


def generate_mock_lcs_region(cds_num_q, cds_num_t, q_domains, t_domains):
    """generate a mock LCS region"""
    query_cds = []
    for i in range(cds_num_q):
        cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
        cds.strand = 1
        cds.orf_num = i
        query_cds.append(cds)

    target_cds = []
    for i in range(cds_num_t):
        cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
        cds.strand = 1
        cds.orf_num = i
        target_cds.append(cds)

    query_domains = []
    for cds, domains in q_domains.items():
        hsps = []
        for domain in domains:
            hsp = bs_hmmer.HSP(query_cds[cds], domain, 100.0, 0, 100)
            hsps.append(hsp)
            query_domains.append(hsp)
        query_cds[cds].hsps = hsps

    target_domains = []
    for cds, domains in t_domains.items():
        hsps = []
        for domain in domains:
            hsp = bs_hmmer.HSP(target_cds[cds], domain, 100.0, 0, 100)
            hsps.append(hsp)
            target_domains.append(hsp)
        target_cds[cds].hsps = hsps

    return query_domains, target_domains


class TestExtendUtilities(unittest.TestCase):
    """Tests for extension utilities"""

    def test_reset_region(self):
        """Test for reset method on a region"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_region(cds_a)
        record_b = generate_mock_region(cds_b)

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        record_pair.comparable_region = bs_comp.ComparableRegion(
            5, 6, 5, 6, 0, 0, 0, 0, False
        )

        # reset the extension
        bs_comp.extend.reset(record_pair)

        expected_comparable_region = bs_comp.ComparableRegion(
            0,
            len(cds_a),
            0,
            len(cds_b),
            0,
            len(list(record_a.get_hsps())),
            0,
            len(list(record_b.get_hsps())),
            False,
        )

        self.assertEqual(expected_comparable_region, record_pair.comparable_region)

    def test_reset_protocluster(self):
        """Test for reset method on a protocluster"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        record_pair.comparable_region = bs_comp.ComparableRegion(
            5, 6, 5, 6, 0, 0, 0, 0, False
        )

        # reset the extension
        bs_comp.extend.reset(record_pair)

        expected_comparable_region = bs_comp.ComparableRegion(
            0,
            len(cds_a),
            0,
            len(cds_b),
            0,
            len(list(record_a.get_hsps())),
            0,
            len(list(record_b.get_hsps())),
            False,
        )

        self.assertEqual(expected_comparable_region, record_pair.comparable_region)

    def test_check_too_short(self):
        """Test for check method when a region is too short and contains no
        biosynthetic genes
        """

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        record_pair.comparable_region = bs_comp.ComparableRegion(
            5, 6, 5, 6, 5, 6, 5, 6, False
        )

        expected_result = (
            False,
            False,
        )
        actual_result = (
            bs_comp.extend.len_check(record_pair, 0.3),
            bs_comp.extend.biosynthetic_check(record_pair),
        )

        self.assertEqual(expected_result, actual_result)

    def test_check_too_short_biosynth(self):
        """Test for check method when the region contains a biosynthetic cds"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        cds_a[5].gene_kind = "biosynthetic"

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        record_pair.comparable_region = bs_comp.ComparableRegion(
            5, 6, 5, 6, 5, 6, 5, 6, False
        )

        expected_result = (
            False,
            True,
        )
        actual_result = (
            bs_comp.extend.len_check(record_pair, 0.3),
            bs_comp.extend.biosynthetic_check(record_pair),
        )

        self.assertEqual(expected_result, actual_result)

    def test_check_pass(self):
        """Test for check method when the region is long enough without a biosynthetic
        cds
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        record_pair.comparable_region = bs_comp.ComparableRegion(
            5, 9, 5, 9, 5, 9, 5, 9, False
        )

        expected_result = (
            True,
            False,
        )
        actual_result = (
            bs_comp.extend.len_check(record_pair, 0.3),
            bs_comp.extend.biosynthetic_check(record_pair),
        )

        self.assertEqual(expected_result, actual_result)


class TestScoreExtend(unittest.TestCase):
    """Tests for score extension"""

    def test_get_target_indexes(self):
        """Tests for getting target domain (cds) indexes"""

        # create mock cds list
        target = []
        for i in range(5):
            cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
            cds.strand = 1
            cds.orf_num = i
            target.append(cds)

        # flip one cds
        target[4].strand = -1

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

        target_domains = []
        for cds in target:
            if cds.strand == 1:
                target_domains.extend(cds.hsps)
            else:
                target_domains.extend(cds.hsps[::-1])

        expected_target_index = {
            "A": [(0, 0), (4, 6)],
            "B": [(1, 1), (3, 3), (3, 4)],
            "C": [(2, 2), (4, 5)],
        }

        actual_target_index = bs_comp.extend.get_target_indexes(target_domains)

        self.assertEqual(expected_target_index, actual_target_index)

    def test_get_target_indexes_reverse(self):
        """Check strand orientation on target with reverse"""
        # create mock cds list
        target = []
        for i in range(5):
            cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
            cds.orf_num = i
            cds.strand = 1
            target.append(cds)

            # flip one cds
        target[4].strand = -1

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

        target_domains = []
        for cds in target:
            if cds.strand == 1:
                target_domains.extend(cds.hsps)
            else:
                target_domains.extend(cds.hsps[::-1])

        # reverse target
        target = target[::-1]
        target_domains = target_domains[::-1]

        expected_index = {
            "A": [(0, 0), (4, 6)],
            "C": [(0, 1), (2, 4)],
            "B": [(1, 2), (1, 3), (3, 5)],
        }

        actual_index = bs_comp.extend.get_target_indexes(target_domains)

        self.assertEqual(expected_index, actual_index)

    def test_get_query_indexes(self):
        """Tests for getting query domain indexes"""
        query = []
        for i in range(3):
            cds = bs_genbank.CDS(i * 100, (i + 1) * 100)
            cds.orf_num = i
            cds.strand = 1
            query.append(cds)

        query[0].hsps = [bs_hmmer.HSP(query[0], "A", 100.0, 0, 100)]
        query[1].hsps = [
            bs_hmmer.HSP(query[1], "B", 100.0, 0, 100),
            bs_hmmer.HSP(query[1], "B", 100.0, 0, 100),
            bs_hmmer.HSP(query[1], "B", 100.0, 0, 100),
        ]
        query[2].hsps = [bs_hmmer.HSP(query[2], "C", 100.0, 0, 100)]

        expected_index = {0: 0, 1: 1, 2: 1, 3: 1, 4: 2}

        query_domains = []
        for cds in query:
            query_domains.extend(cds.hsps)

        actual_index = bs_comp.extend.get_query_indexes(query_domains)

        self.assertEqual(expected_index, actual_index)

    def test_extend_all_match(self):
        """Tests for extension with all matches"""
        # all cds have a common domains
        query, target = generate_mock_cds_lists(
            5, 5, [0, 1, 2, 3, 4], [0, 1, 2, 3, 4], False
        )
        query_domains = []
        for cds in query:
            query_domains.extend(cds.hsps)

        target_domains = []
        for cds in target:
            target_domains.extend(cds.hsps)

        target_index = bs_comp.extend.get_target_indexes(target_domains)
        query_index = bs_comp.extend.get_query_indexes(query_domains)

        expected_extends = (5, 5, 5, 5, 25)

        actual_extends = bs_comp.extend.score_extend(
            query_domains, query_index, 0, 0, target_index, 0, 0, 5, -5, -3, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_large_gap(self):
        """Tests for extension with a large gap"""
        # only first and last have common domains
        # so gap of 3. 3*-2 = -6 with one match = -1. so no extension

        q_domains = {0: ["A"], 1: ["B"]}
        t_domains = {0: ["A"], 1: ["C"], 2: ["D"], 3: ["E"], 4: ["B"]}

        query_dom, target_dom = generate_mock_lcs_region(2, 5, q_domains, t_domains)
        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (1, 1, 1, 1, 5)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 0, 0, target_index, 0, 0, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_small_gap(self):
        """Tests for extension with a small gap"""
        # two gaps, 3 matches. 3 * 5 + 2 * -2 = 11. so extend the full range
        q_domains = {0: ["A"], 1: ["B"], 2: ["C"]}
        t_domains = {0: ["A"], 1: ["B"], 2: ["D"], 3: ["E"], 4: ["C"]}

        query_dom, target_dom = generate_mock_lcs_region(3, 5, q_domains, t_domains)
        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 5, 5, 11)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 0, 0, target_index, 0, 0, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_fill_gap(self):
        """Tests the concept of filling a gap during extend"""
        # the idea is that we have a situation like this:
        # query:  A B C
        # target: A C B
        # when we see something like this, we will have a gap of 1 when we reach B in
        # the query:
        # query:  A - B C
        # target: A C B
        # but we want to "fix" it when we see that close by there is another same domain
        # that would fill the gap

        q_domains = {0: ["A"], 1: ["B"], 2: ["C"]}
        t_domains = {0: ["A"], 1: ["C"], 2: ["B"]}
        query_dom, target_dom = generate_mock_lcs_region(3, 3, q_domains, t_domains)

        # so we expect an extension of 3 on both sides, and a score of 15
        expected_extends = (3, 3, 3, 3, 15)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 0, 0, target_index, 0, 0, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_middle_lcs(self):
        """Tests correct start of extension when lcs is in middle of sequence"""
        q_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        t_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 3, 3, 15)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 2, 2, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_middle_lcs_fill_gap(self):
        """Tests correct gap filling with lcs"""
        q_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        t_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["C"], 4: ["B"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 3, 3, 15)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 2, 2, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_middle_lcs_duplicated_domain(self):
        """Tests correct lcs extension ignores domains before lcs"""
        q_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        t_domains = {0: ["A"], 1: ["X"], 2: ["X"], 3: ["B"], 4: ["C"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        # B and C match, A mismatch -> 5+5-3 = 7
        expected_extends = (3, 3, 2, 2, 7)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 3, 3, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_score_extend_multi_domain_query(self):
        """Tests extension on cds with multiple domains"""
        q_domains = {0: ["X"], 1: ["X"], 2: ["A", "B", "C"]}
        t_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        query_dom, target_dom = generate_mock_lcs_region(3, 5, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (1, 3, 3, 3, 15)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 2, 2, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_score_extend_multi_domain_target(self):
        """Tests extension on cds with multiple domains"""
        q_domains = {0: ["X"], 1: ["X"], 2: ["A"], 3: ["B"], 4: ["C"]}
        t_domains = {0: ["X"], 1: ["X"], 2: ["A", "B", "C"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 3, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 1, 3, 15)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 2, 2, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_score_extend_multidomain_combinations(self):
        """Tests extend on a combination of query and target multidomain cds"""

        q_domains = {1: ["X"], 2: ["X"], 3: ["Q"], 4: ["A", "B"], 5: ["C", "D", "E"]}
        t_domains = {0: ["X"], 1: ["X", "X", "X"], 2: ["A"], 3: ["D", "B"], 4: ["E"]}
        query_dom, target_dom = generate_mock_lcs_region(6, 5, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        # 4 matches ABDE, 2 mismatches QC, : 4*5 - 2*3 = 14
        expected_extends = (3, 6, 3, 4, 14)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 2, target_index, 2, 4, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_forward_mid_cds_target_start(self):
        """Tests extend when starting in the middle of cds"""
        q_domains = {3: ["Q"], 4: ["A", "B"], 5: ["C"]}
        t_domains = {0: ["X"], 1: ["X", "X", "Q"], 2: ["A"], 3: ["B", "C"]}
        query_dom, target_dom = generate_mock_lcs_region(6, 5, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 4, 2, 4, 20)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 0, 0, target_index, 2, 3, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_forward_mid_cds_query_start(self):
        """Tests extend when starting in the middle of cds"""
        q_domains = {2: ["X", "X"], 3: ["X", "Q"], 4: ["A", "B"], 5: ["C"]}
        t_domains = {0: ["X"], 1: ["X", "X", "Q"], 2: ["A"], 3: ["B", "C"]}
        query_dom, target_dom = generate_mock_lcs_region(6, 5, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (2, 4, 2, 4, 20)

        actual_extends = bs_comp.extend.score_extend(
            query_dom, query_index, 2, 3, target_index, 2, 3, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_score_extend_rev(self):
        """Tests correct extension in reverse"""
        q_domains = {0: ["A"], 1: ["B"], 2: ["C"], 3: ["X"], 4: ["X"]}
        t_domains = {0: ["A"], 1: ["B"], 2: ["C"], 3: ["X"], 4: ["X"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 3, 3, 15)

        actual_extends = bs_comp.extend.score_extend_rev(
            query_dom, query_index, 3, 3, target_index, 3, 3, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_rev_no_match_before_lcs(self):
        """Tests correct extension in reverse"""
        q_domains = {0: ["A"], 1: ["B"], 2: ["C"], 3: ["X"], 4: ["X"]}
        t_domains = {0: ["A"], 1: ["B"], 2: ["X"], 3: ["X"], 4: ["C"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        expected_extends = (3, 3, 2, 2, 7)

        actual_extends = bs_comp.extend.score_extend_rev(
            query_dom, query_index, 3, 3, target_index, 2, 2, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_rev_multidomain_combinations(self):
        """Tests extend reverse on a combination of query and target multidomain cds"""
        q_domains = {0: ["Q"], 1: ["E", "A"], 2: ["C", "D"], 3: ["X"], 4: ["X"]}
        t_domains = {0: ["E"], 1: ["D"], 2: ["C", "B"], 3: ["A"], 4: ["X", "X", "X"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        # 4 matches ACDE, 1 gap B: 4*5 - 2 = 18
        expected_extends = (2, 4, 4, 5, 18)

        actual_extends = bs_comp.extend.score_extend_rev(
            query_dom, query_index, 3, 5, target_index, 4, 5, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_rev_double_domains(self):
        """Tests correct extension in reverse with double domains"""
        q_domains = {0: ["A"], 1: ["B", "B"], 2: ["C"], 3: ["X"], 4: ["X"]}
        t_domains = {0: ["A"], 1: ["B"], 2: ["C", "C"], 3: ["X"], 4: ["X"]}
        query_dom, target_dom = generate_mock_lcs_region(5, 6, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        # 3 matches ABC, 1 gap C, 1 mismatch B: 3*5 - 2 - 3 = 10
        expected_extends = (3, 4, 3, 4, 10)

        actual_extends = bs_comp.extend.score_extend_rev(
            query_dom, query_index, 3, 4, target_index, 3, 4, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)

    def test_extend_reverse_mid_cds_starts(self):
        """Tests extend when starting in the middle of cds"""
        q_domains = {2: ["C"], 3: ["B", "A"], 4: ["Q", "X"], 5: ["X", "A"]}
        t_domains = {0: ["B", "C"], 1: ["A"], 2: ["N", "Q", "X"], 3: ["X"]}
        query_dom, target_dom = generate_mock_lcs_region(6, 5, q_domains, t_domains)

        target_index = bs_comp.extend.get_target_indexes(target_dom)
        query_index = bs_comp.extend.get_query_indexes(query_dom)

        # four matches QABC, one gap: 4*5 - 2 = 18
        expected_extends = (2, 4, 2, 5, 18)

        actual_extends = bs_comp.extend.score_extend_rev(
            query_dom, query_index, 2, 4, target_index, 2, 5, 5, -3, -2, 10
        )

        self.assertEqual(expected_extends, actual_extends)


class TestExtendGlocal(unittest.TestCase):
    """Tests for Glocal extension"""

    def test_extend_glocal(self):
        """Tests extend local"""
        # easy case: one short (record A) and one long (record B)
        # should extend both downstream/upstream arms of record A
        #
        # A:           XXXABCXXXXX
        # B: XXXXXXXXXXXXXABCXXXXXXXXX
        #
        a_cds, b_cds = generate_mock_cds_lists(10, 25, [3, 4, 5], [12, 13, 14], False)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)
        pair.comparable_region = bs_comp.ComparableRegion(
            3, 6, 12, 15, 3, 6, 12, 15, False
        )
        bs_comp.extend.extend_glocal(pair)
        expected_glocal = bs_comp.ComparableRegion(0, 10, 12, 15, 0, 10, 12, 15, False)

        conditions = [
            pair.comparable_region == expected_glocal,  # tests cds start/stops
            pair.comparable_region.domain_a_start == expected_glocal.domain_a_start,
            pair.comparable_region.domain_b_start == expected_glocal.domain_b_start,
            pair.comparable_region.domain_a_stop == expected_glocal.domain_a_stop,
            pair.comparable_region.domain_b_stop == expected_glocal.domain_b_stop,
        ]

        self.assertTrue(all(conditions))

    def test_extend_glocal_reverse(self):
        """Tests glocal extend in reverse comparable region"""
        a_cds, b_cds = generate_mock_cds_lists(10, 25, [3, 4, 5], [12, 13, 14], True)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)
        pair.comparable_region = bs_comp.ComparableRegion(
            3, 6, 12, 15, 3, 6, 12, 15, True
        )
        bs_comp.extend.extend_glocal(pair)
        expected_glocal = bs_comp.ComparableRegion(0, 10, 12, 15, 0, 10, 12, 15, True)

        conditions = [
            pair.comparable_region == expected_glocal,  # tests cds start/stops
            pair.comparable_region.domain_a_start == expected_glocal.domain_a_start,
            pair.comparable_region.domain_b_start == expected_glocal.domain_b_start,
            pair.comparable_region.domain_a_stop == expected_glocal.domain_a_stop,
            pair.comparable_region.domain_b_stop == expected_glocal.domain_b_stop,
        ]

        self.assertTrue(all(conditions))

    def test_extend_glocal_diff_arms(self):
        """Tests glocal extend when different record arms are shortest"""
        # both record A and B have a shorter arm
        # should extend upstream A and downstream B arms
        #
        # A:           XXXABCXXXXX
        # B: XXXXXXXXXXXXXABCXX
        #
        a_cds, b_cds = generate_mock_cds_lists(10, 17, [3, 4, 5], [12, 13, 14], False)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)
        pair.comparable_region = bs_comp.ComparableRegion(
            3, 6, 12, 15, 3, 6, 12, 15, False
        )
        bs_comp.extend.extend_glocal(pair)
        expected_glocal = bs_comp.ComparableRegion(0, 6, 12, 17, 0, 6, 12, 17, False)

        conditions = [
            pair.comparable_region == expected_glocal,  # tests cds start/stops
            pair.comparable_region.domain_a_start == expected_glocal.domain_a_start,
            pair.comparable_region.domain_b_start == expected_glocal.domain_b_start,
            pair.comparable_region.domain_a_stop == expected_glocal.domain_a_stop,
            pair.comparable_region.domain_b_stop == expected_glocal.domain_b_stop,
        ]

        self.assertTrue(all(conditions))

    def test_extend_glocal_multi_domain(self):
        """Tests glocal extend with multi domain cdss"""
        # brackets indicate a cds with multiple domains
        #
        # A:      [XX]XX[A BC] DE XXXX
        # B: XXXXXXXXXXXXA[BC][DE]X[XXXX]
        #
        a_cds, b_cds = generate_mock_cds_lists(
            10, 17, [3, 3, 3, 4, 5], [12, 13, 13, 14, 14], False
        )
        a_cds[0].hsps.append(a_cds[0].hsps[0])
        b_cds[-1].hsps.extend([b_cds[-1].hsps[0]] * 3)

        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)
        pair.comparable_region = bs_comp.ComparableRegion(
            3, 6, 12, 15, 4, 9, 12, 17, False
        )
        bs_comp.extend.extend_glocal(pair)
        expected_glocal = bs_comp.ComparableRegion(0, 10, 12, 15, 0, 13, 12, 17, False)

        conditions = [
            pair.comparable_region == expected_glocal,  # tests cds start/stops
            pair.comparable_region.domain_a_start == expected_glocal.domain_a_start,
            pair.comparable_region.domain_b_start == expected_glocal.domain_b_start,
            pair.comparable_region.domain_a_stop == expected_glocal.domain_a_stop,
            pair.comparable_region.domain_b_stop == expected_glocal.domain_b_stop,
        ]

        self.assertTrue(all(conditions))

    def test_extend_greedy(self):
        """Tests greedy extension

        This method of extension should always extend the region to the maximum possible
        CDS that still have a common domain between two records.

        E.g. if we have the following two records:

        A: XAXXBXXXXCX
        B: XXXXXXXAXXXXXBCXXXXXXXXX

        The comparable region should be extended to the following:

        A: XAXXBXXXXCX
            [-------]
        B: XXXXXXXAXXXXXBCXXXXXXXXX
                  [------]
        """

        a_cds, b_cds = generate_mock_cds_lists(11, 24, [1, 4, 9], [11, 13, 14], False)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        bs_comp.extend.extend_greedy(pair)
        expected_greedy = bs_comp.ComparableRegion(1, 10, 11, 15, 1, 10, 11, 15, False)

        conditions = [
            pair.comparable_region == expected_greedy,  # tests cds start/stops
            pair.comparable_region.domain_a_start == expected_greedy.domain_a_start,
            pair.comparable_region.domain_b_start == expected_greedy.domain_b_start,
            pair.comparable_region.domain_a_stop == expected_greedy.domain_a_stop,
            pair.comparable_region.domain_b_stop == expected_greedy.domain_b_stop,
        ]

        self.assertTrue(all(conditions))

    def test_match_extend(self):
        """Tests the new match extend implementation

        This implmentation will be more lenient for domain matches which are
        present on the other side of the lcs. Altogether this should be a more relaxed
        approach to extension and will result in larger comparable regions

        The following example:

        A: XXXXXCXXABXXXX
        B:     XXXXABCXXX
        LCS:       ^^

        Would fail to extend to C in the legacy extend implementation, but should extend
        in the new implementation given a match score of 5 and mismatch of -2
        """

        a_cds, b_cds = generate_mock_cds_lists(14, 10, [6, 9, 10], [5, 6, 7], False)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        # emulate LCS
        pair.comparable_region = bs_comp.ComparableRegion(
            9, 11, 6, 8, 9, 11, 6, 8, False
        )

        bs_comp.extend.extend_simple_match(pair, 5, -2)

        expected_comparable_region = bs_comp.ComparableRegion(
            6, 11, 5, 8, 6, 11, 5, 8, False
        )

        self.assertEqual(pair.comparable_region, expected_comparable_region)
        self.assertEqual(
            pair.comparable_region.domain_a_start,
            expected_comparable_region.domain_a_start,
        )
        self.assertEqual(
            pair.comparable_region.domain_b_start,
            expected_comparable_region.domain_b_start,
        )
        self.assertEqual(
            pair.comparable_region.domain_a_stop,
            expected_comparable_region.domain_a_stop,
        )
        self.assertEqual(
            pair.comparable_region.domain_b_stop,
            expected_comparable_region.domain_b_stop,
        )

    def test_extend_simple_match_multi_domain(self):
        """Tests simple match extend with multi domain cdss"""
        # brackets indicate a cds with multiple domains
        #
        # A:        [XX]BX[A BC]DEX XXXX
        # B: XXXXXXXXXX XX A[BC]EDX[XXXX]
        #
        a_cds, b_cds = generate_mock_cds_lists(
            10, 17, [3, 3, 3, 4, 5], [12, 13, 13, 15, 14], False
        )
        a_cds[0].hsps.append(a_cds[0].hsps[0])
        a_cds[1].hsps = [a_cds[3].hsps[1]]
        b_cds[-1].hsps.extend([b_cds[-1].hsps[0]] * 3)
        record_a = generate_mock_region(a_cds)
        record_b = generate_mock_region(b_cds)
        pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)
        pair.comparable_region = bs_comp.ComparableRegion(
            3, 4, 12, 14, 4, 7, 12, 15, False
        )
        bs_comp.extend.extend_simple_match(pair, 5, -2)
        expected_comparable_region = bs_comp.ComparableRegion(
            1, 6, 12, 16, 2, 9, 12, 17, False
        )
        self.assertEqual(pair.comparable_region, expected_comparable_region)
        self.assertEqual(
            pair.comparable_region.domain_a_start,
            expected_comparable_region.domain_a_start,
        )
        self.assertEqual(
            pair.comparable_region.domain_b_start,
            expected_comparable_region.domain_b_start,
        )
        self.assertEqual(
            pair.comparable_region.domain_a_stop,
            expected_comparable_region.domain_a_stop,
        )
        self.assertEqual(
            pair.comparable_region.domain_b_stop,
            expected_comparable_region.domain_b_stop,
        )
