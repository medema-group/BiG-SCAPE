"""module containing a test for loading files"""

# from python
from unittest import TestCase
from pathlib import Path

# from dependencies
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord

# from other modules
from big_scape.enums import SOURCE_TYPE
import big_scape.enums as bs_enums
from big_scape.genbank import (
    GBK,
    CDS,
    Region,
    ProtoCluster,
    ProtoCore,
    CandidateCluster,
)
from big_scape.cli.config import BigscapeConfig
from big_scape.utility import class_filter, category_filter

# from this module
from big_scape.file_input.load_files import (
    load_dataset_folder,
    load_gbk,
    get_all_bgc_records,
    get_all_bgc_records_query,
    filter_files,
    is_included,
    bgc_length_contraint,
    remove_duplicate_gbk,
    find_mibig_version_url,
)


def create_mock_complete_single_gbk(
    i,
    source_type: bs_enums.SOURCE_TYPE,
    product: list[str] = ["T1PKS"],
    category: list[str] = ["PKS"],
) -> GBK:
    region_feature = SeqFeature(FeatureLocation(0, 100), type="region")
    region_feature.qualifiers = {
        "region_number": ["1"],
        "candidate_cluster_numbers": ["1"],
        "product": product,
    }

    region = Region.parse_as5(region_feature)

    candidate_cluster_feature = SeqFeature(FeatureLocation(0, 100), type="cand_cluster")
    candidate_cluster_feature.qualifiers = {
        "candidate_cluster_number": ["1"],
        "kind": ["single"],
        "protoclusters": ["1"],
        "product": product,
    }

    candidate_cluster = CandidateCluster.parse(candidate_cluster_feature)

    protocluster_feature = SeqFeature(FeatureLocation(0, 100), type="protocluster")
    protocluster_feature.qualifiers = {
        "protocluster_number": ["1"],
        "category": category,
        "product": product,
    }

    protocluster = ProtoCluster.parse(protocluster_feature)

    protocore_feature = SeqFeature(FeatureLocation(0, 100), type="proto_core")
    protocore_feature.qualifiers = {
        "protocluster_number": ["1"],
        "product": product,
    }

    protocore = ProtoCore.parse(protocore_feature)

    protocluster.add_proto_core(protocore)
    candidate_cluster.add_proto_cluster(protocluster)
    region.add_cand_cluster(candidate_cluster)

    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), source_type)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = region
    gbk.region._db_id = i
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    return gbk


def create_mock_gbk(i, sequence, product="test") -> GBK:
    gbk = GBK(Path(f"test_path_{i}.gbk"), str(i), bs_enums.SOURCE_TYPE.QUERY)
    cds = CDS(0, 100)
    cds.parent_gbk = gbk
    cds.orf_num = 1
    cds.strand = 1
    gbk.genes.append(cds)
    gbk.region = Region(gbk, 1, 0, 100, False, product)
    gbk.region._db_id = i
    gbk.metadata = {
        "organism": "banana",
        "taxonomy": "bananus;fruticus",
        "description": "you can eat it",
    }
    seq = SeqRecord(id="test", seq=sequence)
    gbk.nt_seq = seq
    return gbk


class TestLoadGBK(TestCase):
    """Test class for loading of genbank files"""

    def test_gbk_path_invalid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a file, not a folder

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/valid_input.gbk"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
        }

        with self.assertRaises(NotADirectoryError):
            load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

    def test_gbk_path_valid(self):
        """Tests whether loading a given path returns none"""
        # path pointing to a folder containing a valid gbk file

        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "force_gbk": False,
            "include_categories": set(),
            "exclude_categories": set(),
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 4
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_gbk_path_empty(self):
        """Tests whether loading a given path returns none if the folder is empty"""
        # path pointing to a folder containing no gbk files

        run = {
            "input_dir": Path("test/test_data/empty_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        with self.assertRaises(FileNotFoundError):
            load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

    def test_load_gbk_not_a_file(self):
        """Tests whether the load_gbk function correctly returns none when input is not a file"""
        # path pointing to a folder, not a file
        gbk_file_path = Path("test/test_data/valid_gbk_folder/")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        with self.assertRaises(IsADirectoryError):
            load_gbk(gbk_file_path, SOURCE_TYPE.QUERY, run)

    def test_load_gbk_valid(self):
        """Tests whether the load_gbk function correctly loads a valid file

        nb: this does not test the contents of the GBK file, just the existence and reading of it
        """

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        run = {
            "input_dir": Path("test/test_data/valid_gbk_folder/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }
        gbk = load_gbk(gbk_file_path, SOURCE_TYPE.QUERY, run)

        self.assertIsNot(gbk, None)

    def test_include_exclude_str(self):
        """Tests whether include and exclude gbk strings filter correclty"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 2
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_include_all_override_exclude(self):
        """Tests whether include_gbk * correclty includes all files in dir"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": ["*"],
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
        }

        load_result = load_dataset_folder(run["input_dir"], run, SOURCE_TYPE.QUERY)

        expected_count = 4
        actual_count = len(load_result)

        self.assertEqual(expected_count, actual_count)

    def test_filter_files(self):
        """Test test whether a file is included in the include list"""

        files = [
            Path("gbk1_region.gbk"),
            Path("gbk2_region.gbk"),
            Path("gbk3_region.gbk"),
            Path("gbk_region.final.gbk"),
        ]

        include = ["region"]
        exclude = ["final"]

        filtered_files = filter_files(files, include, exclude)

        exp_filtered_files = [
            Path("gbk1_region.gbk"),
            Path("gbk2_region.gbk"),
            Path("gbk3_region.gbk"),
        ]

        self.assertEqual(exp_filtered_files, filtered_files)

    def test_is_included_str(self):
        """Tests whether a filename includes string from list"""

        filename = Path("gbk1_region.gbk")
        include = ["region"]

        self.assertTrue(is_included(filename, include))

    def test_bgc_len_constraint(self):
        """Test whether BGCs are correclty filtered if bellow lenght constraint"""

        gbk_10 = create_mock_gbk(10, "NNNNNNNNNN")
        gbk_20 = create_mock_gbk(20, "NNNNNNNNNNNNNNNNNNNN")
        gbk_30 = create_mock_gbk(30, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")

        gbks = [gbk_10, gbk_20, gbk_30]

        BigscapeConfig.MIN_BGC_LENGTH = 25

        filtered_gbks = bgc_length_contraint(gbks)

        expected_count = 1
        actual_count = len(filtered_gbks)

        self.assertEqual(expected_count, actual_count)

    def test_bgc_len_constraint_error(self):
        """Test whether an error is raised when no BGCs pass the length filter"""

        gbk_10 = create_mock_gbk(10, "NNNNNNNNNN")
        gbk_20 = create_mock_gbk(20, "NNNNNNNNNNNNNNNNNNNN")
        gbk_30 = create_mock_gbk(30, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")

        gbks = [gbk_10, gbk_20, gbk_30]

        BigscapeConfig.MIN_BGC_LENGTH = 35

        self.assertRaises(RuntimeError, bgc_length_contraint, gbks)

    def test_remove_duplicate_gbks(self):
        """Tests wether duplicates are correclty removed"""

        gbk_a = create_mock_gbk(1, "N")
        gbk_b = create_mock_gbk(1, "N")
        gbk_c = create_mock_gbk(2, "NN")

        gbks = [gbk_a, gbk_b, gbk_c]

        filtered_gbks = remove_duplicate_gbk(gbks)

        expected_count = 2
        actual_count = len(filtered_gbks)

        self.assertEqual(expected_count, actual_count)

    def test_include_category(self):
        """Tests whether include category filter works correctly"""
        gbk_a = create_mock_complete_single_gbk(1, SOURCE_TYPE.QUERY, category=["NRPS"])
        gbk_b = create_mock_complete_single_gbk(2, SOURCE_TYPE.QUERY, category=["PKS"])

        # include any records with NRPS
        run = {"include_categories": {"PKS"}, "exclude_categories": set()}

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b]]

        filtered_bgcs = category_filter(run, bgc_records)

        self.assertEqual(filtered_bgcs, [bgc_records[1]])

    def test_exclude_category(self):
        """Tests whether exclude category filter works correctly"""
        gbk_a = create_mock_complete_single_gbk(1, SOURCE_TYPE.QUERY, category=["NRPS"])
        gbk_b = create_mock_complete_single_gbk(
            2, SOURCE_TYPE.QUERY, category=["PKS.NRPS"]  # merged protocluster
        )

        # exclude any records with PKS
        run = {"include_categories": {"NRPS"}, "exclude_categories": {"PKS"}}

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b]]

        filtered_bgcs = category_filter(run, bgc_records)

        self.assertEqual(filtered_bgcs, [bgc_records[0]])

    def test_include_class(self):
        """Tests whether include class filter works correctly"""
        gbk_a = create_mock_gbk(1, "N", product="NRPS.PKS")
        gbk_b = create_mock_gbk(2, "N", product="NRPS")
        gbk_c = create_mock_gbk(3, "N", product="NRPS-like")

        # include any records with NRPS
        run = {"include_classes": {"NRPS"}, "exclude_classes": set()}

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b, gbk_c]]

        filtered_bgcs = class_filter(run, bgc_records)

        # NRPS-like is not included
        self.assertEqual(filtered_bgcs, bgc_records[0:2])

    def test_exclude_class(self):
        """Tests whether exclude class filter works correctly"""
        gbk_a = create_mock_gbk(1, "N", product="NRPS.PKS")
        gbk_b = create_mock_gbk(2, "N", product="NRPS")
        gbk_c = create_mock_gbk(3, "N", product="NRPS-like")

        # exclude any records with NRPS
        run = {"include_classes": set(), "exclude_classes": {"NRPS"}}

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b, gbk_c]]

        filtered_bgcs = class_filter(run, bgc_records)

        # NRPS-like is not excluded
        self.assertEqual(filtered_bgcs, [bgc_records[-1]])

    def test_include_and_exclude_class(self):
        """Tests whether include and exclude class filter work together correctly"""
        gbk_a = create_mock_gbk(1, "N", product="NRPS.PKS")
        gbk_b = create_mock_gbk(2, "N", product="NRPS")
        gbk_c = create_mock_gbk(3, "N", product="NRPS-like")

        # exclude any PKS, then only keep NRPS
        run = {"include_classes": {"NRPS"}, "exclude_classes": {"PKS"}}

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b, gbk_c]]

        filtered_bgcs = class_filter(run, bgc_records)

        # PKS is excluded, meaning only the pure NRPS is kept
        self.assertEqual(filtered_bgcs, [bgc_records[1]])

    def test_category_and_class_filter(self):
        """Tests whether category filter works together with class correctly"""
        gbk_a = create_mock_complete_single_gbk(
            1, SOURCE_TYPE.QUERY, product=["NRPS"], category=["NRPS"]
        )
        gbk_b = create_mock_complete_single_gbk(
            2,
            SOURCE_TYPE.QUERY,
            product=["NRPS", "T1PKS"],
            category=["PKS.NRPS"],  # merged protocluster
        )
        gbk_c = create_mock_complete_single_gbk(
            3,
            SOURCE_TYPE.QUERY,
            product=["NRPS", "T2PKS"],
            category=["PKS.NRPS"],  # merged protocluster
        )

        # include all PKS category except T1PKS products
        run = {
            "include_categories": {"PKS"},
            "exclude_categories": set(),
            "include_classes": set(),
            "exclude_classes": {"T1PKS"},
        }

        bgc_records = [gbk.region for gbk in [gbk_a, gbk_b, gbk_c]]

        categ_bgcs = category_filter(run, bgc_records)
        filtered_bgcs = class_filter(run, categ_bgcs)

        self.assertEqual(filtered_bgcs, [bgc_records[2]])

    def test_get_all_bgc_records(self):
        """Tests whether all records are correctly loaded from a list of gbks"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "record_type": bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER,
        }

        gbks = []

        for i in range(3):
            gbk = create_mock_complete_single_gbk(i, SOURCE_TYPE.QUERY)
            gbks.append(gbk)

        all_bgc_records = get_all_bgc_records(run, gbks)

        expected_count = 3
        actual_count = len(all_bgc_records)

        self.assertEqual(expected_count, actual_count)

    def test_get_all_bgc_records_query(self):
        """Tests whether all records are correctly loaded from a list of gbks, as well as query record"""

        run = {
            "input_dir": Path("test/test_data/alt_valid_gbk_input/"),
            "input_mode": bs_enums.INPUT_MODE.RECURSIVE,
            "include_gbk": None,
            "exclude_gbk": None,
            "cds_overlap_cutoff": None,
            "cores": None,
            "classify": False,
            "record_type": bs_enums.genbank.RECORD_TYPE.PROTO_CLUSTER,
            "query_record_number": 1,
        }

        gbks = []

        for i in range(3):
            gbk = create_mock_complete_single_gbk(i, SOURCE_TYPE.REFERENCE)
            gbks.append(gbk)

        query_gbk = create_mock_complete_single_gbk(3, SOURCE_TYPE.QUERY)
        gbks.append(query_gbk)

        exp_query_record = query_gbk.region.cand_clusters[1].proto_clusters[1]

        all_bgc_records, query_record = get_all_bgc_records_query(run, gbks)

        self.assertEqual(exp_query_record, query_record)

    def test_find_mibig_version_url(self):
        """Tests correct fetching of mibig download urls"""

        self.skipTest("Reduce sending requests")

        version = "3.1"
        url_31 = find_mibig_version_url(version)
        self.assertEqual(
            "https://dl.secondarymetabolites.org/mibig/mibig_antismash_3.1_gbk.tar.bz2",
            url_31,
        )

        version = "4.0"
        url_40 = find_mibig_version_url(version)
        self.assertEqual(
            "https://dl.secondarymetabolites.org/mibig/mibig_antismash_4.0_gbk_as8b1.tar.bz2",
            url_40,
        )

        version = "1"
        self.assertRaises(ValueError, find_mibig_version_url, version)
