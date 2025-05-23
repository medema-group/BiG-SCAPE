"""Contains tests for the GBK component"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.dereplicating.gbk_components.gbk import GBK
import big_scape.enums as bs_enums
from big_scape.dereplicating.gbk_component_parsing import validate_cds_component
from big_scape.cli.config import BigscapeConfig


class TestGBKComponent(TestCase):
    """Test class for the GBK component"""

    def test_create(self):
        """Tests whether the create function correctly creates a GBK object"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk_hash = "hash"
        nt_length = 100
        as_version = "1.0"
        source_type = bs_enums.SOURCE_TYPE.QUERY

        gbk = GBK(gbk_file_path, gbk_hash, nt_length, as_version, source_type)

        self.assertIsInstance(gbk, GBK)

    def test_validate_CDS(self):
        """Tests whether validate cds component correclty returns false if none present"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk_hash = "hash"
        nt_length = 100
        as_version = "1.0"
        source_type = bs_enums.SOURCE_TYPE.QUERY

        gbk = GBK(gbk_file_path, gbk_hash, nt_length, as_version, source_type)

        self.assertFalse(validate_cds_component(gbk))

    def test_bgc_len_filter_below(self):
        """Tests whether BGCs are correctly filtered out by length constraint"""

        gbk_10 = GBK(Path("gbk_path"), "gbk_hash", 10, "as_version", bs_enums.SOURCE_TYPE.QUERY)
        gbk_20 = GBK(Path("gbk_path"), "gbk_hash", 20, "as_version", bs_enums.SOURCE_TYPE.QUERY)
        gbk_30 = GBK(Path("gbk_path"), "gbk_hash", 30, "as_version", bs_enums.SOURCE_TYPE.QUERY)

        gbk_list = [gbk_10, gbk_20, gbk_30]

        BigscapeConfig.MIN_BGC_LENGTH = 25

        filtered_gbks = GBK.length_filter(gbk_list)

        self.assertEqual(len(filtered_gbks), 1)

    def test_bgc_len_filter_none_left(self):
        """Tests whether BGCs are correctly filtered out by length constraint"""

        gbk_10 = GBK(Path("gbk_path"), "gbk_hash", 10, "as_version", bs_enums.SOURCE_TYPE.QUERY)
        gbk_20 = GBK(Path("gbk_path"), "gbk_hash", 20, "as_version", bs_enums.SOURCE_TYPE.QUERY)
        gbk_30 = GBK(Path("gbk_path"), "gbk_hash", 30, "as_version", bs_enums.SOURCE_TYPE.QUERY)

        gbk_list = [gbk_10, gbk_20, gbk_30]

        BigscapeConfig.MIN_BGC_LENGTH = 35

        self.assertRaises(RuntimeError, GBK.length_filter, gbk_list)

