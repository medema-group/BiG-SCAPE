"""Contains tests for the lcs module."""

# from python
import unittest
import random
import string

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.hmm as bs_hmmer
import big_scape.comparison as bs_comparison


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


class TestCDSLCS(unittest.TestCase):
    """Tests for the lcs module using CDS lists."""

    def test_lcs_cds_forward(self):
        """Test lcs detection for two CDS lists"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 2, 3], [1, 2, 3], False)

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 1, 4, False))

    def test_lcs_cds_reverse(self):
        """Test lcs detection for two CDS lists, reverse"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 2, 3], [1, 2, 3], True)

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 6, 9, True))

    def test_lcs_cds_empty(self):
        """Test lcs detection for two empty CDS lists"""
        cds_a, cds_b = generate_mock_cds_lists(0, 0, [], [], False)

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (0, 0, 0, 0, False))

    def test_lcs_cds_len_one(self):
        """Test lcs detection for two CDS lists where there are only matches of len=1"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], False)

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 2, 1, 2, False))

    def test_lcs_cds_len_one_reverse(self):
        """Test lcs detection for two CDS lists where there are only matches of len=1,
        reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], True)

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 2, 8, 9, False))

    def test_lcs_cds_len_multiple(self):
        """Test lcs detection for two CDS lists where there are equal matches of len>1"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 3, 5, 6], [1, 2, 3, 5, 6], False
        )

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 1, 4, False))

    def test_lcs_cds_len_multiple_reverse(self):
        """Test lcs detection for two CDS lists where there are equal matches of len>1,
        reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 3, 5, 6], [1, 2, 3, 5, 6], True
        )

        lcs = bs_comparison.lcs.find_cds_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 4, 6, 9, True))


class testDomainLCS(unittest.TestCase):
    """Tests for the lcs module using domain lists."""

    def test_lcs_domains_forward(self):
        """Test lcs detection for two domain lists"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], False
        )

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 6, 1, 6, False))

    def test_lcs_domains_reverse(self):
        """Test lcs detection for two domain lists, reverse"""
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 2, 3], [1, 2, 2, 2, 3], True
        )

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 6, 6, 11, True))

    def test_lcs_domains_empty(self):
        """Test lcs detection for two empty domain lists"""
        cds_a, cds_b = generate_mock_cds_lists(0, 0, [], [], False)

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (0, 0, 0, 0, False))

    def test_lcs_domains_len_one(self):
        """Test lcs detection for two domain lists where there are only matches of len=1"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 5], False)

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        # should return most central, not first
        self.assertEqual(lcs, (5, 6, 5, 6, False))

    def test_lcs_domains_len_one_reverse(self):
        """Test lcs detection for two domain lists where there are only matches of
        len=1, reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [1, 3, 5], [1, 3, 6], True)

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        # should return most central, not first
        self.assertEqual(lcs, (1, 2, 8, 9, False))

    def test_lcs_domains_len_multiple(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], False
        )

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 5, 1, 5, False))

    def test_lcs_domains_len_multiple_reverse(self):
        """Test lcs detection for two domain lists where there are equal matches of
        len>1, reverse
        """
        cds_a, cds_b = generate_mock_cds_lists(
            10, 10, [1, 2, 2, 3, 5, 6], [1, 2, 2, 3, 5, 6], True
        )

        lcs = bs_comparison.lcs.find_domain_lcs(cds_a, cds_b)

        self.assertEqual(lcs, (1, 5, 6, 10, True))
