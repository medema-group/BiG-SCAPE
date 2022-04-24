"""Module containing bigslice filter parameter helper functions and classes

Author: Arjan Draisma
"""

from src import utility


class BigsliceParam():
    """Class which keeps track of run options relating to bigslice pre-filter"""
    # whether to use mibig in a run
    use_bigslice = False

    # where to store any files downloaded for the functionality
    bigslice_data_folder: str

    # cutoff value to use for distances
    bigslice_cutoff: float

    ANTISMASH_URL = "https://github.com/" + \
        "antismash/antismash/archive/5-1-1.tar.gz"
    ANTISMASH_VERSION = "antismash-5-1-1"
    REFERENCE_PROTEINS_URL = "ftp://ftp.pir.georgetown.edu/" + \
        "databases/rps/rp-seqs-15.fasta.gz"

    # build_subpfam algorithm parameters
    MIN_CLADE_NUMBERS = 3
    MAX_CLADE_NUMBERS = 50
    MIN_PROTEIN_SEQUENCES = 5
    
    def __init__(self, options):
        self.set_bigslice(options)
        self.set_bigslice_folder(options)
        return

    def set_bigslice(self, options):
        """Set basic values relating to bigslice"""
        self.use_bigslice = options.bigslice_filter
        self.bigslice_cutoff = options.bs_filter_cutoff
        return

    def set_bigslice_folder(self, options):
        """Sets the bigslice data folder and creates it if it does not yet exist"""
        self.bigslice_data_folder = options.bigslice_data_folder
        utility.create_directory(self.bigslice_data_folder, "bigslice_data", False)
        return
