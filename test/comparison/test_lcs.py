"""Contains tests for the lcs module."""

# from python
import unittest
import random
import string

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.hmm as bs_hmmer
import big_scape.comparison as bs_comparison
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


def generate_mock_protocluster(cds_list: list[bs_genbank.CDS], protocore_idx: int):
    """Generate a mock protocluster from a cds list"""
    gbk = bs_genbank.GBK(None, bs_enums.SOURCE_TYPE.QUERY)
    gbk.genes = cds_list
    protocluster = bs_genbank.ProtoCluster(
        gbk, 1, 0, len(cds_list) * 100, False, "", "", {}
    )

    protocore_start = cds_list[protocore_idx].nt_start
    protocore_stop = cds_list[protocore_idx].nt_stop

    protocore = bs_genbank.ProtoCore(gbk, 1, protocore_start, protocore_stop, False, "")

    protocluster.proto_core[1] = None
    protocluster.add_proto_core(protocore)

    return protocluster


class TestCDSLCS(unittest.TestCase):
    """Tests for the lcs module using CDS lists."""

    def test_lcs_cds_forward(self):
        """Test lcs detection for two CDS lists"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 2, 3], [1, 2, 3], False)

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 1, 4, False))

    def test_lcs_cds_reverse(self):
        """Test lcs detection for two CDS lists, reverse"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 2, 3], [1, 2, 3], True)

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 6, 9, True))

    def test_lcs_cds_empty(self):
        """Test lcs detection for two empty CDS lists"""
        cds_a, cds_b = generate_mock_cds_lists(0, 0, [], [], False)

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (0, 0, 0, 0, False))

    def test_lcs_cds_len_one(self):
        """Test lcs detection for two CDS lists where there are only matches of len=1"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], False)

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 2, 1, 2, False))

    def test_lcs_cds_len_one_reverse(self):
        """Test lcs detection for two CDS lists where there are only matches of len=1,
        reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], True)

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 2, 8, 9, False))

    def test_lcs_cds_len_multiple(self):
        """Test lcs detection for two CDS lists where there are equal matches of len>1"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 3, 5, 6], [1, 2, 3, 5, 6], False
        )

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 1, 4, False))

    def test_lcs_cds_len_multiple_reverse(self):
        """Test lcs detection for two CDS lists where there are equal matches of len>1,
        reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 3, 5, 6], [1, 2, 3, 5, 6], True
        )

        lcs = bs_comparison.lcs.find_cds_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 6, 9, True))


class testDomainLCS(unittest.TestCase):
    """Tests for the lcs module using domain lists."""

    def test_lcs_domains_forward(self):
        """Test lcs detection for two domain lists"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 6, 1, 6, False))

    def test_lcs_domains_reverse(self):
        """Test lcs detection for two domain lists, reverse"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], True
        )

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 6, 6, 11, True))

    def test_lcs_domains_empty(self):
        """Test lcs detection for two empty domain lists"""
        cds_a, cds_b = generate_mock_cds_lists(0, 0, [], [], False)

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (0, 0, 0, 0, False))

    def test_lcs_domains_len_one(self):
        """Test lcs detection for two domain lists where there are only matches of len=1"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], False)

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        # should return most central, not first
        self.assertEqual(lcs, (5, 6, 5, 6, False))

    def test_lcs_domains_len_one_reverse(self):
        """Test lcs detection for two domain lists where there are only matches of
        len=1, reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 6], True)

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        # should return most central, not first
        self.assertEqual(lcs, (1, 2, 8, 9, False))

    def test_lcs_domains_len_multiple(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], False
        )

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 5, 1, 5, False))

    def test_lcs_domains_len_multiple_reverse(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1, reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], True
        )

        lcs = bs_comparison.lcs.find_domain_lcs_region(cds_a, cds_b)

        self.assertEqual(lcs, (1, 5, 6, 10, True))


class TestProtoclusterDomainLCS(unittest.TestCase):
    """Tests for the lcs module using protocluster lists."""

    def test_find_protocore_distance_before(self):
        """Tests whether the distance to a protocore from a given index is calculated
        correctly
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        # protocore is third cds
        protocluster = generate_mock_protocluster(cds_a, 2)

        expected_distance = 2

        actual_distance = bs_comparison.lcs.find_protocore_distance(protocluster, 0)

        self.assertEqual(expected_distance, actual_distance)

    def test_find_protocore_distance_in_protocore(self):
        """Tests whether the distance to a protocore from a given index is calculated
        correctly
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        # protocore is third cds
        protocluster = generate_mock_protocluster(cds_a, 2)

        expected_distance = 0

        actual_distance = bs_comparison.lcs.find_protocore_distance(protocluster, 2)

        self.assertEqual(expected_distance, actual_distance)

    def test_find_protocore_distance_after(self):
        """Tests whether the distance to a protocore from a given index is calculated
        correctly
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        # protocore is third cds
        protocluster = generate_mock_protocluster(cds_a, 2)

        expected_distance = 5

        actual_distance = bs_comparison.lcs.find_protocore_distance(protocluster, 7)

        self.assertEqual(expected_distance, actual_distance)
