#!/usr/bin/env python


"""
BiG-SCAPE

PI: Marnix Medema               marnix.medema@wur.nl

Maintainers/developers:
Jorge Navarro                   j.navarro@westerdijkinstitute.nl
Satria Kautsar                  satria.kautsar@wur.nl

Developer:
Emmanuel (Emzo) de los Santos   E.De-Los-Santos@warwick.ac.uk


Usage:   Please see `python bigscape.py -h`

Example: python bigscape.py -c 8 --pfam_dir ./ -i ./inputfiles -o ./results

Status: beta

Official repository:
https://git.wageningenur.nl/medema-group/BiG-SCAPE


# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.
"""

# Makes sure the script can be used with Python 2 as well as Python 3.
from __future__ import print_function
from __future__ import division


import os
import logging
import warnings
import time

from distutils import dir_util

from sys import version_info

# refactored imports
from src import big_scape
from src import gbk
from src import hmm
from src import mibig
from src import pfam
from src import utility

# disable pyling complaining about backwards compatibility measures
#pylint: disable=wrong-import-order,redefined-builtin,invalid-name,undefined-variable,import-error
if version_info[0] == 2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0] == 3:
    import pickle # for storing and retrieving dictionaries
#pylint: enable=wrong-import-order,redefined-builtin,invalid-name,undefined-variable,import-error

def init_logger(root_path, options):
    """Initializes the logger for big-scape

    input:
        root_path - root path of this bigscape.py file
        options - options object containing command line arument information
    """

    ## logging
    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")
    root_logger = logging.getLogger()

    # create log dir
    log_dir = os.path.join(root_path, "log")
    if not os.path.exists(log_dir):
        os.mkdir(log_dir)
    # set log file
    log_time_stamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    log_file = "log/" + log_time_stamp + ".log"

    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)

    if not options.quiet:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(log_formatter)
        root_logger.addHandler(console_handler)

    if options.verbose:
        root_logger.level = logging.DEBUG
    else:
        root_logger.level = logging.INFO

if __name__ == "__main__":
    # get root path of this project
    ROOT_PATH = os.path.dirname(__file__)

    # get run options
    # ROOT_PATH is passed here because the imports no longer allow us to use __file__
    OPTIONS = utility.cmd_parser(ROOT_PATH)

    # init logger
    init_logger(ROOT_PATH, OPTIONS)

    # ignore specific warnings from sklearn
    # TODO: investigate whether this can be mediated by changing the affinity propagation.
    # concerns the following warning:
    # UserWarning: All samples have mutually equal similarities. Returning arbitrary cluster
    # center(s).
    warnings.filterwarnings(action="ignore", category=UserWarning)

    # create new run details
    # ideally we parse all the options and remember them in this object
    # we will later use options to describe what parameters were used, and
    # use the run object to describe what effect that had, e.g. which files
    # were used
    RUN = big_scape.Run()

    # this should be the only occasion where options is needed. no further than this.
    RUN.init(OPTIONS)

    # preparation is done. run timing formally starts here
    RUN.start()


    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    logging.info("   - - Processing input files - -")

    mibig.extract_mibig(RUN)

    # there are three steps where BiG-SCAPE will import GBK files:
    # 1. mibig
    # 2. genbank
    # 3. query BGC
    # all of these are generalized into this method, which returns the bgc info and gen bank info
    BGC_INFO_DICT, GBK_FILE_DICT, MIBIG_SET = gbk.import_gbks(RUN)
    
    # base name as a list and as a set
    # also used in all vs all analysis
    CLUSTER_NAME_LIST = list(GBK_FILE_DICT.keys())
    CLUSTER_NAME_SET = set(CLUSTER_NAME_LIST)


    ### Step 2: Run hmmscan
    logging.info("Trying threading on %d cores", RUN.options.cores)
    logging.info("Predicting domains using hmmscan")

    # get all fasta files in cache directory
    CACHED_FASTA_FILES = hmm.get_cached_fasta_files(RUN)

    # verify that all clusters have a corresponding fasta file in cache
    hmm.check_fasta_files(RUN, CLUSTER_NAME_SET, CACHED_FASTA_FILES)

    # Make a list of all fasta files that need to be processed by hmmscan & BiG-SCAPE
    # (i.e., they don't yet have a corresponding .domtable)
    # includes all fasta files if force_hmmscan is set
    FASTA_FILES_TO_PROCESS = hmm.get_fasta_files_to_process(RUN, CACHED_FASTA_FILES)

    # if any are there, run hmmscan
    if len(FASTA_FILES_TO_PROCESS) > 0:
        # this function blocks the main thread until finished
        hmm.run_hmmscan_async(RUN, FASTA_FILES_TO_PROCESS)
        logging.info(" Finished generating domtable files.")
    else:
        logging.info(" All files were processed by hmmscan. Skipping step...")


    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    logging.info("Parsing hmmscan domtable files")

    # All available domtable files
    CACHED_DOMTABLE_FILES = hmm.get_cached_domtable_files(RUN)

    # verify that domtable files were generated successfully. each cluster should have a domtable
    # file.
    hmm.check_domtable_files(RUN, CLUSTER_NAME_SET, CACHED_DOMTABLE_FILES)

    # find unprocessed files (assuming that if the pfd file exists, the pfs should too)
    # this will just return all domtable files if force_hmmscan is set
    DOMTABLE_FILES_TO_PROCESS = hmm.get_domtable_files_to_process(RUN, CACHED_DOMTABLE_FILES)

    for domtableFile in DOMTABLE_FILES_TO_PROCESS:
        hmm.parse_hmmscan(domtableFile, RUN.directories.pfd, RUN.directories.pfs,
                          RUN.options.domain_overlap_cutoff, RUN.options.verbose, GBK_FILE_DICT,
                          CLUSTER_NAME_LIST, CLUSTER_NAME_SET, MIBIG_SET)

    # If number of pfd files did not change, no new sequences were added to the
    #  domain fastas and we could try to resume the multiple alignment phase
    # baseNames have been pruned of BGCs with no domains that might've been added temporarily
    TRY_RESUME_MULTIPLE_ALIGNMENT = False
    ALREADY_DONE = hmm.get_processed_domtable_files(RUN, CACHED_DOMTABLE_FILES)
    if len(CLUSTER_NAME_SET - set(pfd.split(os.sep)[-1][:-9] for pfd in ALREADY_DONE)) == 0:
        TRY_RESUME_MULTIPLE_ALIGNMENT = True
    else:
        # new sequences will be added to the domain fasta files. Clean domains folder
        # We could try to make it so it's not necessary to re-calculate every alignment,
        #  either by expanding previous alignment files or at the very least,
        #  re-aligning only the domain files of the newly added BGCs
        logging.info(" New domain sequences to be added; cleaning domains folder")
        for thing in os.listdir(RUN.directories.domains):
            os.remove(os.path.join(RUN.directories.domains, thing))

    logging.info(" Finished generating pfs and pfd files.")


    ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per
    ### sample objects
    logging.info("Processing domains sequence files")

    # do one more check of pfd files to see if they are all there
    hmm.check_pfd_files(RUN, CLUSTER_NAME_SET)

    # BGCs --
    # this collection will contain all bgc objects
    BGC_COLLECTION = big_scape.BgcCollection()

    # init the collection with the acquired names from importing GBK files
    BGC_COLLECTION.initialize(CLUSTER_NAME_LIST)

    # the BGCs need to know which domains belong where
    # this is done in this step
    BGC_COLLECTION.init_ordered_domain_list(RUN)

    # add info and source gbk files
    BGC_COLLECTION.add_bgc_info(BGC_INFO_DICT)
    BGC_COLLECTION.add_source_gbk_files(GBK_FILE_DICT)

    if RUN.options.skip_ma:
        logging.info(" Running with skip_ma parameter: Assuming that the domains folder has all the \
            fasta files")
        BGC_COLLECTION.load_domain_names_from_dict(RUN)
    else:
        logging.info(" Adding sequences to corresponding domains file")
        BGC_COLLECTION.load_domain_names_from_pfd(RUN, TRY_RESUME_MULTIPLE_ALIGNMENT)

        BGC_COLLECTION.save_domain_names_to_dict(RUN)

    # Key: BGC. Item: ordered list of simple integers with the number of domains
    # in each gene
    # Instead of `DomainCountGene = defaultdict(list)`, let's try arrays of
    # unsigned ints
    # GENE_DOMAIN_COUNT = {}
    # list of gene-numbers that have a hit in the anchor domain list. Zero based
    # COREBIOSYNTHETIC_POS = {}
    # list of +/- orientation
    # BGC_GENE_ORIENTATION = {}

    # TODO: remove this comment? not sure what it relates to
    # if it's a re-run, the pfd/pfs files were not changed, so the skip_ma flag
    # is activated. We have to open the pfd files to get the gene labels for
    # each domain
    # We now always have to have this data so the alignments are produced
    PARSE_PFD_RESULTS = big_scape.parse_pfd(RUN, BGC_COLLECTION)
    GENE_DOMAIN_COUNT, COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION = PARSE_PFD_RESULTS

    BGC_COLLECTION.add_gene_domain_counts(GENE_DOMAIN_COUNT)
    BGC_COLLECTION.add_bio_synth_core_pos(COREBIOSYNTHETIC_POS)
    BGC_COLLECTION.add_gene_orientations(BGC_GENE_ORIENTATION)

    # at this point we can assemble a gene string for bgc info
    BGC_COLLECTION.init_gene_strings()



    ### Step 5: Create SVG figures
    logging.info(" Creating arrower-like figures for each BGC")

    # read hmm file. We'll need that info anyway for final visualization
    logging.info("  Parsing hmm file for domain information")
    PFAM_INFO = pfam.parse_pfam_a(RUN)
    logging.info("    Done")

    # verify if there are figures already generated

    big_scape.generate_images(RUN, CLUSTER_NAME_SET, GBK_FILE_DICT, PFAM_INFO, BGC_INFO_DICT)
    logging.info(" Finished creating figures")
    logging.info("   - - Calculating distance matrix - -")

    # Do multiple alignments if needed
    if not RUN.options.skip_ma:
        hmm.do_multiple_align(RUN, TRY_RESUME_MULTIPLE_ALIGNMENT)

    # If there's something to analyze, load the aligned sequences
    logging.info(" Trying to read domain alignments (*.algn files)")
    ALIGNED_DOMAIN_SEQS = hmm.read_aligned_files(RUN)

    # copy html templates
    HTML_TEMPLATE_PATH = os.path.join(ROOT_PATH, "html_template", "output")
    dir_util.copy_tree(HTML_TEMPLATE_PATH, RUN.directories.output)

    # make a new run folder in the html output & copy the overview_html
    big_scape.copy_template_per_cutoff(RUN, ROOT_PATH)
    RUNDATA_NETWORKS_PER_RUN = big_scape.prepare_cutoff_rundata_networks(RUN)
    HTML_SUBS_PER_RUN = big_scape.prepare_html_subs_per_run(RUN)

    # create pfams.js
    pfam.create_pfam_js(RUN, PFAM_INFO)

    # Try to make default analysis using all files found inside the input folder
    logging.info("Generating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    logging.info("   Writing the complete Annotations file for the complete set")
    big_scape.write_network_annotation_file(RUN, BGC_COLLECTION)

    # Find index of all MIBiG BGCs if necessary
    if RUN.mibig.use_mibig:
        NAME_TO_IDX = {}
        for clusterIdx, clusterName in enumerate(BGC_COLLECTION.bgc_name_tuple):
            NAME_TO_IDX[clusterName] = clusterIdx

        MIBIG_SET_INDICES = set()
        for bgc in MIBIG_SET:
            MIBIG_SET_INDICES.add(NAME_TO_IDX[bgc])

    # Making network files mixing all classes
    if RUN.options.mix:
        big_scape.generate_network(RUN, BGC_COLLECTION, ALIGNED_DOMAIN_SEQS,
                                   MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                   HTML_SUBS_PER_RUN, True)

    # Making network files separating by BGC class
    if not RUN.options.no_classify:
        big_scape.generate_network(RUN, BGC_COLLECTION, ALIGNED_DOMAIN_SEQS,
                                   MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                   HTML_SUBS_PER_RUN, False)

    # fetch genome list for overview.js
    INPUT_CLUSTERS_IDX = []
    big_scape.fetch_genome_list(RUN, INPUT_CLUSTERS_IDX, BGC_COLLECTION.bgc_name_tuple, MIBIG_SET, BGC_INFO_DICT,
                                GBK_FILE_DICT)

    # update family data (convert global bgc indexes into input-only indexes)
    big_scape.update_family_data(RUNDATA_NETWORKS_PER_RUN, INPUT_CLUSTERS_IDX, BGC_COLLECTION.bgc_name_tuple,
                                 MIBIG_SET)

    # generate overview data
    RUN.end()
    big_scape.generate_results_per_cutoff_value(RUN, RUNDATA_NETWORKS_PER_RUN, HTML_SUBS_PER_RUN)
    # dump bgc info
    pickle.dump(BGC_INFO_DICT, open(os.path.join(RUN.directories.cache, 'bgc_info.dict'), 'wb'))

    # done
    RUN.report_runtime()
