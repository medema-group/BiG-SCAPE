"""Contains tests for the lcs module."""

# from python
import unittest
import random
import string
import big_scape.comparison.record_pair

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.hmm as bs_hmmer
import big_scape.comparison as bs_comparison
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
            cds_a.append(a)
        if i < len_b:
            b = bs_genbank.CDS(i * 100, (i + 1) * 100)
            b.strand = 1
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


def generate_mock_protocluster(cds_list: list[bs_genbank.CDS], protocore_idx: int = -1):
    """Generate a mock protocluster from a cds list"""
    gbk = bs_genbank.GBK(None, "test", bs_enums.SOURCE_TYPE.QUERY)
    gbk.genes = cds_list
    protocluster = bs_genbank.ProtoCluster(
        gbk, 1, 0, len(cds_list) * 100, False, "", {}
    )

    if protocore_idx > -1:
        protocore_start = cds_list[protocore_idx].nt_start
        protocore_stop = cds_list[protocore_idx].nt_stop

        protocore = bs_genbank.ProtoCore(
            gbk, 1, protocore_start, protocore_stop, False, ""
        )

        protocluster.proto_core[1] = None
        protocluster.add_proto_core(protocore)

    return protocluster


class TestRegionDomainLCS(unittest.TestCase):
    """Tests for the lcs module using region domain lists."""

    def test_lcs_domains_forward(self):
        """Test lcs detection for two domain lists

        should return the longest match, not the most central
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (1, 6, 1, 6, 1, 4, 1, 4, False))

    def test_lcs_domains_reverse(self):
        """Test lcs detection for two domain lists, reverse

        should return the longest match, not the most central
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], True
        )

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (1, 6, 1, 6, 1, 4, 1, 4, True))

    def test_lcs_domains_empty(self):
        """Test lcs detection for two empty domain lists

        this should throw a runtime error"""
        cds_a, cds_b = generate_mock_cds_lists(0, 0, [], [], False)

        pc_a = generate_mock_protocluster(cds_a)
        pc_b = generate_mock_protocluster(cds_b)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        with self.assertRaises(RuntimeError):
            bs_comparison.lcs.find_domain_lcs_region(pair)

    def test_lcs_domains_len_one(self):
        """Test lcs detection for two domain lists where there are only matches of len=1

        should return the most central match, not the first
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], False)

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        # should return most central, not first
        self.assertEqual(lcs, (5, 6, 5, 6, 5, 6, 5, 6, False))

    def test_lcs_domains_len_one_reverse(self):
        """Test lcs detection for two domain lists where there are only matches of
        len=1, reverse

        should return the most central match, not the first
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], True)

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (5, 6, 5, 6, 5, 6, 5, 6, True))

    def test_lcs_domains_len_multiple(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], False
        )

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (1, 5, 1, 5, 1, 4, 1, 4, False))

    def test_lcs_domains_len_multiple_reverse(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1, reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], True
        )

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (1, 5, 1, 5, 1, 4, 1, 4, True))

    def test_lcs_domains_len_multiple_unaligned(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [2, 3, 3, 4, 6, 7], [1, 2, 2, 3, 5, 6], False
        )

        pc_a = generate_mock_protocluster(cds_a, 2)
        pc_b = generate_mock_protocluster(cds_b, 2)
        pair = big_scape.comparison.record_pair.RecordPair(pc_a, pc_b)

        lcs = bs_comparison.lcs.find_domain_lcs_region(pair)

        self.assertEqual(lcs, (2, 6, 1, 5, 2, 5, 1, 4, False))


class TestProtoclusterDomainLCS(unittest.TestCase):
    """Tests for the lcs module using protocluster lists."""

    def test_lcs_protocluster_type_check_a(self):
        """Tests whether the correct type error is raised when record_a in a recordpair
        is not a protocluster
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        gbk = bs_genbank.GBK("", "", bs_enums.SOURCE_TYPE.QUERY)

        record_a = bs_genbank.Region(gbk, 1, 0, len(cds_a) * 100, False, "")
        record_b = bs_genbank.Region(gbk, 1, 0, len(cds_b) * 100, False, "")

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        with self.assertRaises(TypeError):
            bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

    def test_lcs_protocluster_type_check_b(self):
        """Tests whether the correct type error is raised when record_b in a recordpair
        is not a protocluster
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        gbk = bs_genbank.GBK("", "", bs_enums.SOURCE_TYPE.QUERY)

        record_a = bs_genbank.ProtoCluster(
            gbk, 1, 0, len(cds_a) * 100, False, "", "", {}
        )
        record_b = bs_genbank.Region(gbk, 1, 0, len(cds_b) * 100, False, "")

        record_pair = big_scape.comparison.record_pair.RecordPair(record_a, record_b)

        with self.assertRaises(TypeError):
            bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

    def test_lcs_protocluster_forward(self):
        """Test lcs detection for two domain lists"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        protocluster_a = generate_mock_protocluster(cds_a, 2)
        protocluster_b = generate_mock_protocluster(cds_b, 2)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (1, 6, 1, 6, 1, 4, 1, 4, False)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_reverse(self):
        """Test lcs detection for two domain lists, reverse"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], True
        )

        protocluster_a = generate_mock_protocluster(cds_a, 2)
        protocluster_b = generate_mock_protocluster(cds_b, 2)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (1, 6, 1, 6, 1, 4, 1, 4, True)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_eq_len_fwd_protocore(self):
        """Test lcs detection for two protoclusters. Equal length matches forward and
        reverse, but one of the matches has both cds in protocore

        In this case the forward match is in the protocore, so should be returned
        """

        # cds lists:
        #     |-->|   |-->|
        # a: 0 1 2 3 4 5 6 7 8 9
        # b: 0 1 2 3 4 5 6 7 8 9
        #     |-->|   |-->|
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 5, 6], [1, 2, 5, 6], False
        )

        # cds lists:
        #     |-->|   |-->|
        # a: 0 1 2 3 4 5 6 7 8 9
        # b: 0 1 2 3 4 5 6 7 8 9
        #     |-->|   |<--|
        cds_b[5:7] = cds_b[5:7][::-1]

        protocluster_a = generate_mock_protocluster(cds_a, 2)
        protocluster_b = generate_mock_protocluster(cds_b, 2)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (1, 3, 1, 3, 1, 3, 1, 3, False)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_eq_len_rev_protocore(self):
        """Test lcs detection for two protoclusters. Equal length matches forward and
        reverse, but one of the matches has both cds in protocore

        In this case the reverse match is in the protocore and should be returned
        """

        # cds lists:
        #     |-->|   |-->|
        # a: 0 1 2 3 4 5 6 7 8 9
        # b: 0 1 2 3 4 5 6 7 8 9
        #     |-->|   |-->|
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 5, 6], [1, 2, 5, 6], False
        )

        # cds lists:
        #     |-->|   |-->|
        # a: 0 1 2 3 4 5 6 7 8 9
        # b: 0 1 2 3 4 5 6 7 8 9
        #     |-->|   |<--|
        cds_b[5:7] = cds_b[5:7][::-1]

        protocluster_a = generate_mock_protocluster(cds_a, 5)
        protocluster_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (5, 7, 3, 5, 5, 7, 3, 5, True)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_eq_len_1_protocore(self):
        """Test lcs detection for two protoclusters. Single-length matches but one of
        the matches has both cds in protocore
        """

        # cds lists:
        # PC          |-|
        #     |>|     |>|   |>|
        # a: 0 1 2 3 4 5 6 7 8 9
        # b: 0 1 2 3 4 5 6 7 8 9
        #     |>|     |>| |>|
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 5, 8], [1, 5, 7], False)

        protocluster_a = generate_mock_protocluster(cds_a, 5)
        protocluster_b = generate_mock_protocluster(cds_b, 5)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (5, 6, 5, 6, 5, 6, 5, 6, False)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_rev_fwd_protocore_fwd_longer(self):
        """Test lcs detection for two protoclusters. This is a case where a match is
        found of the same length forward and reverse, and they're both in the protocore

        This prefers the forward match unless the reverse match is longer. in this case,
        the forward is longer
        """

        # domain/cds lists:
        # CDS: 0 1 2 3 4         5 6
        #  PC:        |---------|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |-------->|
        cds_a, cds_b = generate_mock_cds_lists(
            7, 7, [1, 4, 4, 4, 4, 4], [1, 4, 4, 4, 4, 4], False
        )

        # domain/cds lists:
        # CDS: 0 1 2 3 4         5 6
        #  PC:        |---------|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |---->|<--|
        cds_b[4].hsps[3:5] = cds_b[4].hsps[3:5][::-1]

        protocluster_a = generate_mock_protocluster(cds_a, 4)
        protocluster_b = generate_mock_protocluster(cds_b, 4)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (4, 7, 4, 7, 4, 5, 4, 5, False)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_rev_fwd_protocore_fwd_same(self):
        """Test lcs detection for two protoclusters. This is a case where a match is
        found of the same length forward and reverse, and they're both in the protocore

        This prefers the match closer to the center of the shortest domainlist. Since the lists are
        of the same length, this defaults to A. In any case, the forward match should be chosen.
        """

        # domain/cds lists:
        # CDS: 0 1 2 3 4       5 6 7
        #  PC:        |-------|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |------>|
        cds_a, cds_b = generate_mock_cds_lists(
            8, 8, [1, 4, 4, 4, 4], [1, 4, 4, 4, 4], False
        )

        # domain/cds lists:
        # CDS: 0 1 2 3 4       5 6 7
        #  PC:        |-------|
        #       |>|   |------>|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |-->|<--|
        cds_b[4].hsps[2:4] = cds_b[4].hsps[2:4][::-1]

        protocluster_a = generate_mock_protocluster(cds_a, 4)
        protocluster_b = generate_mock_protocluster(cds_b, 4)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (4, 6, 4, 6, 4, 5, 4, 5, False)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)

    def test_lcs_protocluster_rev_fwd_protocore_rev_longer(self):
        """Test lcs detection for two protoclusters. This is a case where there are two
        matches in the protocore, but the reverse is longer and should be preferred
        """

        # domain/cds lists:
        # CDS: 0 1 2 3 4         5 6
        #  PC:        |---------|
        #       |>|   |-------->|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |-------->|
        cds_a, cds_b = generate_mock_cds_lists(
            7, 7, [1, 4, 4, 4, 4, 4], [1, 4, 4, 4, 4, 4], False
        )

        # domain/cds lists:
        # CDS: 0 1 2 3 4         5 6
        #  PC:        |---------|
        #       |>|   |-------->|
        #   a: 0 1 2 3 4 5 6 7 8 9 10
        #   b: 0 1 2 3 4 5 6 7 8 9 10
        #       |>|   |-->|<----|
        cds_b[4].hsps[2:5] = cds_b[4].hsps[2:5][::-1]

        protocluster_a = generate_mock_protocluster(cds_a, 4)
        protocluster_b = generate_mock_protocluster(cds_b, 4)

        record_pair = big_scape.comparison.record_pair.RecordPair(
            protocluster_a, protocluster_b
        )

        expected_lcs = (6, 9, 2, 5, 4, 5, 2, 3, True)
        actual_lcs = bs_comparison.lcs.find_domain_lcs_protocluster(record_pair)

        self.assertEqual(expected_lcs, actual_lcs)
