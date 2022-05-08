#!/usr/bin/env python3


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
import pickle

from distutils import dir_util

from sys import version_info

# refactored imports
from src import big_scape
from src import gbk
from src import hmm
from src import mibig
from src import pfam
from src import utility
from src import data

def init_logger(options):
    """Initializes the logger for big-scape

    input:
        root_path - root path of this bigscape.py file
        options - options object containing command line arument information
    """

    ## logging
    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")
    root_logger = logging.getLogger()

    # create log dir
    if not os.path.exists(options.log_path):
        os.mkdir(options.log_path)
    # set log file
    log_time_stamp = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
    log_file = os.path.join(options.log_path, log_time_stamp + ".log")

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
    ## Initialization steps
    # show an error if we're not on python 3
    if version_info.major != 3:
        logging.error("This python program was developed for Python 3. Please make sure to use \
            python 3 when running BiG-SCAPE")
        exit(1)

    # get root path of this project
    ROOT_PATH = os.path.realpath(os.path.dirname(__file__))

    # get run options
    # ROOT_PATH is passed here because the imports no longer allow us to use __file__
    OPTIONS = utility.cmd_parser(ROOT_PATH)

    # init logger
    init_logger(OPTIONS)

    # initialize the profiler & start it
    PROFILER = utility.Profiler(OPTIONS)
    PROFILER.start()

    # ignore specific warnings from sklearn
    # TODO: investigate whether this can be mediated by changing the affinity propagation.
    # concerns the following warning:
    # UserWarning: All samples have mutually equal similarities. Returning arbitrary cluster
    # center(s).
    warnings.filterwarnings(action="ignore", category=UserWarning)

    # TODO: download whatever version of pfam is necessary if it isn't there

    # TODO: next press the hmm file into more optimal formats
    # this currently does not work since it generates hmm files that don't work with the
    # optimized profiles used later on. User still needs to generate their own pressed hmms
    # hmm.pyhmmer_hmmpress(OPTIONS)

    # create new run details
    # ideally we parse all the options and remember them in this object
    # we will later use options to describe what parameters were used, and
    # use the run object to describe what effect that had, e.g. which files
    # were used
    RUN = big_scape.Run()

    # this should be the only occasion where options is needed. no further than this.
    RUN.init(OPTIONS)

    logging.info("Options set to use %d cores", RUN.options.cores)

    # preparation is done. run timing formally starts here
    RUN.start()

    # init database
    DB_PATH = os.path.join(RUN.directories.output, "data.db")
    DB = data.Database(DB_PATH, True)

    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    logging.info("   - - Processing input files - -")

    if RUN.mibig.use_mibig:
        mibig.extract_mibig(RUN)

    # load the input data into the database
    data.initialize_db(RUN, DB)

    # base name as a list and as a set
    # also used in all vs all analysis
    BGC_IDS = data.get_cluster_id_list(DB)
    logging.info("Working with %d total bgcs", len(BGC_IDS))

    ### Step 2: Run hmmscan
    # find which sequences have already had their domains predicted
    PREDICTED_IDS = data.get_predicted_bgc_list(DB)

    # find the difference
    IDS_TODO = list(set(BGC_IDS) - set(PREDICTED_IDS))
    # IDS_TODO = list()

    # run hmmscan if there are any sequences which have not yet been processed by hmmscan
    if len(IDS_TODO) > 0:
        logging.info("Predicting %d domains using hmmsearch", len(IDS_TODO))
        # this function blocks the main thread until finished
        hmm.run_pyhmmer_pfam(RUN, DB, IDS_TODO)
        logging.info(" Finished predicting domains")
    else:
        logging.info(" All files were processed by hmmscan. Skipping step...")

    # dumping db here since this is the largest task until this point.
    DB.dump_db_file()

    if RUN.options.feature_filter:
        logging.info(" Generating features")

        features = data.Features.extract(data.get_cluster_id_list(DB), DB)

        for feature in features:
            feature.save(DB)

        DB.commit_inserts()

    # get a list of high scoring protein hits
    HSPS = data.get_hsp_id_list(DB)

    # get a list of already-aligned hsps
    # TODO: work this into a status table like with bgcs
    # currently this just checks if there is an entry for this HSP in the hsp_alignment table
    # TODO: identify changes in dataset and realign
    ALIGNED_HSPS = data.get_aligned_hsp_list(DB)

    # get the alignments to be done
    HSPS_TODO = list(set(HSPS) - set(ALIGNED_HSPS))
    # HSPS_TODO = list()

    # if there are any to be done, we'll align
    if len(HSPS_TODO) > 0:
        logging.info("Performing multiple alignments using hmmalign")
        hmm.do_multiple_align(RUN, DB, HSPS)
    else:
        logging.info(" All high scoring protein domains were already aligned. Skipping step...")


    DB.dump_db_file()

    ### Step 4: Create SVG figures
    logging.info(" Creating arrower-like figures for each BGC")

    # read hmm file. We'll need that info anyway for final visualization
    logging.info("  Parsing hmm file for domain information")
    PFAM_INFO = pfam.parse_pfam_a(RUN)
    logging.info("    Done")

    # verify if there are figures already generated

    big_scape.generate_images(RUN, DB, PFAM_INFO)
    logging.info(" Finished creating figures")
    logging.info("   - - Calculating distance matrix - -")

    # copy html templates
    HTML_TEMPLATE_PATH = os.path.join(ROOT_PATH, "html_template", "output")
    dir_util.copy_tree(HTML_TEMPLATE_PATH, RUN.directories.output)

    # make a new run folder in the html output & copy the overview_html
    big_scape.copy_template_per_cutoff(RUN, ROOT_PATH)
    RUNDATA_NETWORKS_PER_RUN = big_scape.prepare_cutoff_rundata_networks(RUN)
    HTML_SUBS_PER_RUN = big_scape.prepare_html_subs_per_run(RUN)

    # TODO: refactor ArrowerSVG such that it can use the database so we don't have
    # to generate the PFS files. currently arrowersvg requires these files to work
    data.generate_pfd_files(RUN, DB)

    # create pfams.js
    pfam.create_pfam_js(RUN, PFAM_INFO)

    # Try to make default analysis using all files found inside the input folder
    logging.info("Generating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    logging.info("   Writing the complete Annotations file for the complete set")
    # big_scape.write_network_annotation_file(RUN, BGC_COLLECTION)

    # TODO: rework data storage from this point onwards. Now we're converting
    # back to large amounts of data in memory because refactoring the storage
    # after this point is a massive task
    BGC_INFO_DICT, GBK_FILE_DICT, MIBIG_SET = gbk.import_gbks(RUN)
    BGC_COLLECTION = data.generate_bgc_collection(RUN, DB, BGC_INFO_DICT, GBK_FILE_DICT)
    ALIGNED_DOMAIN_SEQS = data.generate_aligned_domain_seqs(RUN, DB)
    MIBIG_SET_INDICES = data.generate_mibig_set_indices(RUN, BGC_COLLECTION, MIBIG_SET)

    # Making network files mixing all classes
    if RUN.options.mix:
        big_scape.generate_network(RUN, DB, BGC_COLLECTION, ALIGNED_DOMAIN_SEQS,
                                   MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                   HTML_SUBS_PER_RUN, True)

    # Making network files separating by BGC class
    if not RUN.options.no_classify:
        big_scape.generate_network(RUN, DB, BGC_COLLECTION, ALIGNED_DOMAIN_SEQS,
                                   MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                   HTML_SUBS_PER_RUN, False)

    # fetch genome list for overview.js
    INPUT_CLUSTERS_IDX = []
    BGC_INFO_DICT = data.gen_bgc_info_for_fetch_genome(DB)
    GBK_FILE_DICT = data.get_cluster_gbk_dict(RUN, DB)
    big_scape.fetch_genome_list(RUN, INPUT_CLUSTERS_IDX, BGC_COLLECTION.bgc_name_tuple, MIBIG_SET,
                                BGC_INFO_DICT, GBK_FILE_DICT)

    # update family data (convert global bgc indexes into input-only indexes)
    big_scape.update_family_data(RUNDATA_NETWORKS_PER_RUN, INPUT_CLUSTERS_IDX,
                                 BGC_COLLECTION.bgc_name_tuple, MIBIG_SET)

    # generate overview data
    RUN.end()
    big_scape.generate_results_per_cutoff_value(RUN, RUNDATA_NETWORKS_PER_RUN, HTML_SUBS_PER_RUN)
    # dump bgc info
    pickle.dump(BGC_INFO_DICT, open(os.path.join(RUN.directories.cache, 'bgc_info.dict'), 'wb'))

    # done
    RUN.report_runtime()

    PROFILER.stop()
