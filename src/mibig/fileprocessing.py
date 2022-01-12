"""Class to support file processing for MIBIG files

Author: Arjan Draisma
"""

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

    print("\n Trying to read bundled MIBiG BGCs as reference")
    print("Assuming mibig path: {}".format(run.options.mibig_path))

    # try to see if the zip file has already been decompressed
    numbgcs = len(glob(os.path.join(run.mibig.gbk_path, "*.gbk")))
    if numbgcs == 0:
        if not zipfile.is_zipfile(run.mibig.zip_path):
            sys.exit("Did not find file {}. \
                Please re-download it from the official repository".format(run.mibig.zip_path))

        with zipfile.ZipFile(run.mibig.zip_path, 'r') as mibig_zip:
            for fname in mibig_zip.namelist():
                if fname[-3:] != "gbk":
                    continue

                extractedbgc = mibig_zip.extract(fname, path=run.options.mibig_path)
                if run.options.verbose:
                    print("  Extracted {}".format(extractedbgc))

    elif run.mibig.expected_num_bgc == numbgcs:
        print("  MIBiG BGCs seem to have been extracted already")
    else:
        sys.exit("Did not find the correct number of MIBiG BGCs ({}). \
            Please clean the 'Annotated MIBiG reference' folder from any \
            .gbk files first".format(run.mibig.expected_num_bgc))


def import_mibig_gbk(run: big_scape.Run, bgc_info, genbank_dict):
    """Imports MIBiG GBK files and stores information into dedicated objects

    Inputs:
    - run: parameters relevant to the current run
    - bgc_info: Stores, per BGC: predicted type, gbk Description, number of records, width of
      longest record, GenBank's accession, Biosynthetic Genes' ids
    - genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    """
    # set of mibig bgcs
    mibig_set = set()

    print("\nImporting MIBiG files")
    mibig_gbk_path = run.mibig.gbk_path
    gbk.get_gbk_files(mibig_gbk_path, run, bgc_info, genbank_dict, True)
    # gbk.get_gbk_files(run.mibig.gbk_path, run.directories.output, run.directories.bgc_fasta,
    #                   int(options.min_bgc_size), ['*'], run.gbk.exclude, bgc_info, options.mode,
    #                   options.verbose, options.force_hmmscan, run.valid_classes,
    #                   bgctools.BgcData, genbank_dict)

    for key in genbank_dict:
        mibig_set.add(key)

    return mibig_set
