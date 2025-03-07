"""Contains tests for the GBK component"""

# from python
from unittest import TestCase
from pathlib import Path

# from other modules
from big_scape.dereplicating.gbk_components.gbk import GBK
import big_scape.enums as bs_enums
from big_scape.dereplicating.gbk_component_parsing import validate_cds_component


class TestGBKComponent(TestCase):
    """Test class for the GBK component"""

    def test_create(self):
        """Tests whether the create function correctly creates a GBK object"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk_hash = "hash"
        nt_length = 100
        as_version = "1.0"
        source_type = bs_enums.SOURCE_TYPE.QUERY

        gbk = GBK.create(gbk_file_path, gbk_hash, nt_length, as_version, source_type)

        self.assertIsInstance(gbk, GBK)

    def test_validate_CDS(self):
        """Tests whether validate cds component correclty returns false if none present"""

        gbk_file_path = Path("test/test_data/valid_gbk_folder/valid_input_region.gbk")
        gbk_hash = "hash"
        nt_length = 100
        as_version = "1.0"
        source_type = bs_enums.SOURCE_TYPE.QUERY

        gbk = GBK.create(gbk_file_path, gbk_hash, nt_length, as_version, source_type)

        self.assertFalse(validate_cds_component(gbk))
