"""Class to support file processing for MIBIG files

Author: Arjan Draisma
"""

import os
import sys
import zipfile
from glob import glob

import src.gbk as gbk
import src.bgctools as bgctools

from src.big_scape.run.base import Run



def read_mibig(run: Run, options):
    """Reads MiBIG files and generates a set of BGCs


    Inputs:
    - mibig_param: parameters relevant to mibig functionality of this run
    - verbose: whether to report code execution verbosely
    """
    # Read included MIBiG
    # Change this for every officially curated MIBiG bundle
    # (file, final folder, number of bgcs)

    
    # Stores, per BGC: predicted type, gbk Description, number of records, width of longest record,
    # GenBank's accession, Biosynthetic Genes' ids
    bgc_info = {}
    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    genbankDict = {}
    # set of mibig bgcs
    mibig_set = set()


    print("\n Trying to read bundled MIBiG BGCs as reference")
    #TODO: __file__ usage, simplify?
    
    mibig_zip_file = run.mibig.mibig_file + ".zip"
    print("Assuming mibig path: {}".format(options.mibig_path))
    bgcs_path = os.path.join(options.mibig_path, run.mibig.mibig_file)

    # try to see if the zip file has already been decompressed
    numbgcs = len(glob(os.path.join(bgcs_path, "*.gbk")))
    if numbgcs == 0:
        if not zipfile.is_zipfile(os.path.join(options.mibig_path, mibig_zip_file)):
            sys.exit("Did not find file {}. \
                Please re-download it from the official repository".format(mibig_zip_file))

        with zipfile.ZipFile(os.path.join(options.mibig_path, mibig_zip_file), 'r') as mibig_zip:
            for fname in mibig_zip.namelist():
                if fname[-3:] != "gbk":
                    continue

                extractedbgc = mibig_zip.extract(fname, path=options.mibig_path)
                if options.verbose:
                    print("  Extracted {}".format(extractedbgc))

    elif run.mibig.expected_num_bgc == numbgcs:
        print("  MIBiG BGCs seem to have been extracted already")
    else:
        sys.exit("Did not find the correct number of MIBiG BGCs ({}). \
            Please clean the 'Annotated MIBiG reference' folder from any \
            .gbk files first".format(run.mibig.expected_num_bgc))

    print("\nImporting MIBiG files")
    gbk.get_gbk_files(bgcs_path, run.directories.output, run.directories.bgc_fasta,
                      int(options.min_bgc_size), ['*'], run.gbk.exclude, bgc_info, options.mode,
                      options.verbose, options.force_hmmscan, run.valid_classes,
                      bgctools.bgc_data, genbankDict)

    for i in genbankDict.keys():
        mibig_set.add(i)

    return mibig_set, bgc_info, genbankDict
