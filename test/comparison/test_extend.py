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
    gbk = bs_genbank.GBK(None, bs_enums.SOURCE_TYPE.QUERY)
    gbk.genes = cds_list

    region = bs_genbank.Region(gbk, 1, 0, len(cds_list) * 100, False, "")

    return region


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


class TestExtend(unittest.TestCase):
    """Tests for extension of comparableregions with region pairs"""

    def test_reset_expansion_region(self):
        """Test for reset_expansion method on a region"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_region(cds_a)
        record_b = generate_mock_region(cds_b)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        # reset the expansion
        bs_comp.extend.reset_expansion(comparable_region)

        expected_comparable_region = bs_comp.ComparableRegion(
            record_pair, 0, len(cds_a), 0, len(cds_b), False
        )

        self.assertEqual(expected_comparable_region, comparable_region)

    def test_reset_expansion_protocluster(self):
        """Test for reset_expansion method on a protocluster"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        # reset the expansion
        bs_comp.extend.reset_expansion(comparable_region)

        expected_comparable_region = bs_comp.ComparableRegion(
            record_pair, 0, len(cds_a), 0, len(cds_b), False
        )

        self.assertEqual(expected_comparable_region, comparable_region)

    def test_check_expand_too_short(self):
        """Test for check_expand method"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        expected_result = False
        actual_result = bs_comp.extend.check_expand(comparable_region)

        self.assertEqual(expected_result, actual_result)

    def test_check_expand_too_short_biosynth(self):
        """Test for check_expand method when the region contains a biosynthetic
        domain"""

        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        cds_a[5].gene_kind = "biosynthetic"

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 6, 5, 6, False)

        expected_result = True
        actual_result = bs_comp.extend.check_expand(comparable_region)

        self.assertEqual(expected_result, actual_result)

    def test_check_expand_pass(self):
        """Test for check_expand method when the region is long enough without
        a biosynthetic cds"""
        cds_a, cds_b = generate_mock_cds_lists(10, 10, [0, 1, 2], [0, 1, 2], False)

        record_a = generate_mock_protocluster(cds_a, 5)
        record_b = generate_mock_protocluster(cds_b, 5)

        record_pair = bs_comp.RecordPair(record_a, record_b)

        # create a comparable region starting at a single cds
        comparable_region = bs_comp.ComparableRegion(record_pair, 5, 9, 5, 9, False)

        expected_result = True
        actual_result = bs_comp.extend.check_expand(comparable_region)

        self.assertEqual(expected_result, actual_result)

    def test_check_lcs_too_short(self):
        """Test for check_lcs method when lcs is too short"""
        self.fail("Not implemented")

    def test_check_lcs_too_short_biosynth(self):
        """Test for check_lcs method when lcs is too short, but it contiains a
        biosynthetic domain"""
        self.fail("Not implemented")

    def test_check_lcs_pass(self):
        """Test for check_lcs method when lcs is long enough"""
        self.fail("Not implemented")

    def test_expand(self):
        """Test for expand method"""
        self.fail("Not implemented")
