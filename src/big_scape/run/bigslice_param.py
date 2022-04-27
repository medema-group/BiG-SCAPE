"""Module containing bigslice filter parameter helper functions and classes

Author: Arjan Draisma
"""

from os import path
from src import utility


class BigsliceParam():
    """Class which keeps track of run options relating to bigslice pre-filter"""
    # whether to use mibig in a run
    use_bigslice = False

    # where to store any files downloaded for the functionality
    bigslice_data_path: str

    # cutoff value to use for distances
    bigslice_cutoff: float

    bigslice_db_location = "https://github.com/medema-group/bigslice/releases/download/v1.0.0/bigslice-models.2020-04-27.tar.gz"
    bigslice_db_md5 = "d875dd0166f1c79167c53970314c3ac7"

    # url to download antismash files from
    antismash_url = "https://github.com/" + \
    "antismash/antismash/archive/5-1-1.tar.gz"
    # this needs to reflect the version number in the tar because of how the tarfile is structured
    antismash_tar_folder = "antismash-5-1-1"

    
    def __init__(self, options):
        self.set_bigslice(options)
        self.set_bigslice_folder(options)
        return

    def set_bigslice(self, options):
        """Set basic values relating to bigslice"""
        self.use_bigslice = options.bigslice_filter
        self.bigslice_cutoff = options.bigslice_filter_cutoff
        return

    def set_bigslice_folder(self, options):
        """Sets the bigslice data folder and creates it if it does not yet exist"""
        self.bigslice_data_path = options.bigslice_data_path
        utility.create_directory(self.bigslice_data_path, "bigslice_data", False)

        biosyn_pfam_dir = path.join(
            self.bigslice_data_path,
            "biosynthetic_pfams"
        )
        utility.create_directory(biosyn_pfam_dir, "bigslice_data", False)
        return
