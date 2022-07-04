"""Module containing functionality to  to support file processing for MIBIG
files

Author: Jorge Navarro, Arjan Draisma
"""

import logging
import os
import sys
import zipfile
from glob import glob

import src.gbk as gbk
import src.big_scape as big_scape

def extract_mibig(run: big_scape.Run):
    """Extracts MIBiG zips

    Inputs:
    - run: parameters relevant to the current run
    - verbose: whether to report code execution verbosely
    """
    # Read included MIBiG
    # TODO: automatically download new versions

    logging.info(" Trying to read bundled MIBiG BGCs as reference")
    logging.info("  Assuming mibig path: %s", run.options.mibig_path)

    # try to see if the zip file has already been decompressed
    numbgcs = len(glob(os.path.join(run.mibig.gbk_path, "*.gbk")))
    if numbgcs == 0:
        if not zipfile.is_zipfile(run.mibig.zip_path):
            logging.error("Did not find file %s. \
                          Please re-download it from the official repository", run.mibig.zip_path)
            sys.exit(1)

        with zipfile.ZipFile(run.mibig.zip_path, 'r') as mibig_zip:
            for fname in mibig_zip.namelist():
                if fname[-3:] != "gbk":
                    continue

                extractedbgc = mibig_zip.extract(fname, path=run.options.mibig_path)

                logging.debug("  Extracted %s", extractedbgc)

    elif run.mibig.expected_num_bgc == numbgcs:
        logging.info("  MIBiG BGCs seem to have been extracted already")
    else:
        logging.error("Did not find the correct number of MIBiG BGCs (%d). \
            Please clean the 'Annotated MIBiG reference' folder from any \
            .gbk files first", run.mibig.expected_num_bgc)
        sys.exit(1)
