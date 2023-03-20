"""Contains methods to press an input .hmm file into more optimal formats for hmmscan and hmmalign
"""

# from python
from pathlib import Path
import logging

# from dependencies
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmpress


def hmm_press(hmm_path: Path):
    """Presses hmm files for optimized files for further analysis

    Raises a ValueError if the given HMM file cannot be read by pyhmmer"""
    with HMMFile(hmm_path) as hmm_file:
        if hmm_file.is_pressed():
            logging.info("PFAM Hmm file was already pressed.")
            return

        logging.info("Pressing Pfam HMM file...")
        hmmpress(hmm_file, hmm_path)
