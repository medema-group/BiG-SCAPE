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
import time

from distutils import dir_util

# refactored imports
from src import big_scape
from src import gbk
from src import hmm
from src import mibig
from src import pfam
from src import utility

from sys import version_info
if version_info[0] == 2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0] == 3:
    import pickle # for storing and retrieving dictionaries

if __name__ == "__main__":
    # get root path of this project
    ROOT_PATH = os.path.basename(__file__)

    # get run options
    # ROOT_PATH is passed here because the imports no longer allow us to use __file__
    OPTIONS = utility.CMD_parser(ROOT_PATH)

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
    print("\n\n   - - Processing input files - -")

    mibig.extract_mibig(RUN)

    # there are three steps where BiG-SCAPE will import GBK files:
    # 1. mibig
    # 2. genbank
    # 3. query BGC
    # but all of them end up calling get_gbk_files, which add onto two dicts:
    # BGC_INFO and GEN_BANK_DICT
    # TODO: generalize into object
    BGC_INFO = {}
    GEN_BANK_DICT = {}

    MIBIG_SET = mibig.import_mibig_gbk(RUN, BGC_INFO, GEN_BANK_DICT)

    print("\nImporting GenBank files")
    gbk.fileprocessing.import_genbank_gbk(RUN, BGC_INFO, GEN_BANK_DICT)

    QUERY_BGC = ""
    if RUN.directories.has_query_bgc:
        print("\nImporting query BGC files")
        QUERY_BGC = gbk.fileprocessing.import_query_gbk(RUN, BGC_INFO, GEN_BANK_DICT)

    # CLUSTERS and SAMPLE_DICT contain the necessary structure for all-vs-all and sample analysis
    CLUSTERS = list(GEN_BANK_DICT.keys())

    print("\nTrying threading on {} cores".format(str(RUN.options.cores)))


    ### Step 2: Run hmmscan
    print("\nPredicting domains using hmmscan")

    CLUSTER_BASE_NAMES = set(CLUSTERS)

    # get all fasta files in cache directory
    CACHED_FASTA_FILES = hmm.get_cached_fasta_files(RUN)

    # verify that all clusters have a corresponding fasta file in cache
    hmm.check_fasta_files(RUN, CLUSTER_BASE_NAMES, CACHED_FASTA_FILES)

    # Make a list of all fasta files that need to be processed by hmmscan & BiG-SCAPE
    # (i.e., they don't yet have a corresponding .domtable)
    # includes all fasta files if force_hmmscan is set
    FASTA_FILES_TO_PROCESS = hmm.get_fasta_files_to_process(RUN, CACHED_FASTA_FILES)

    # if any are there, run hmmscan
    if len(FASTA_FILES_TO_PROCESS) > 0:
        # this function blocks the main thread until finished
        hmm.run_hmmscan_async(RUN, FASTA_FILES_TO_PROCESS)
        print(" Finished generating domtable files.")
    else:
        print(" All files were processed by hmmscan. Skipping step...")


    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    print("\nParsing hmmscan domtable files")

    # All available domtable files
    CACHED_DOMTABLE_FILES = hmm.get_cached_domtable_files(RUN)

    # verify that domtable files were generated successfully. each cluster should have a domtable
    # file.
    hmm.check_domtable_files(RUN, CLUSTER_BASE_NAMES, CACHED_DOMTABLE_FILES)

    # find unprocessed files (assuming that if the pfd file exists, the pfs should too)
    # this will just return all domtable files if force_hmmscan is set
    DOMTABLE_FILES_TO_PROCESS = hmm.get_domtable_files_to_process(RUN, CACHED_DOMTABLE_FILES)

    for domtableFile in DOMTABLE_FILES_TO_PROCESS:
        hmm.parse_hmmscan(domtableFile, RUN.directories.pfd, RUN.directories.pfs,
                          RUN.options.domain_overlap_cutoff, RUN.options.verbose, GEN_BANK_DICT,
                          CLUSTERS, CLUSTER_BASE_NAMES, MIBIG_SET)

    # If number of pfd files did not change, no new sequences were added to the
    #  domain fastas and we could try to resume the multiple alignment phase
    # baseNames have been pruned of BGCs with no domains that might've been added temporarily
    TRY_RESUME_MULTIPLE_ALIGNMENT = False
    ALREADY_DONE = hmm.get_processed_domtable_files(RUN, CACHED_DOMTABLE_FILES)
    if len(CLUSTER_BASE_NAMES - set(pfd.split(os.sep)[-1][:-9] for pfd in ALREADY_DONE)) == 0:
        TRY_RESUME_MULTIPLE_ALIGNMENT = True
    else:
        # new sequences will be added to the domain fasta files. Clean domains folder
        # We could try to make it so it's not necessary to re-calculate every alignment,
        #  either by expanding previous alignment files or at the very least,
        #  re-aligning only the domain files of the newly added BGCs
        print(" New domain sequences to be added; cleaning domains folder")
        for thing in os.listdir(RUN.directories.domains):
            os.remove(os.path.join(RUN.directories.domains, thing))

    print(" Finished generating pfs and pfd files.")


    ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per
    ### sample objects
    print("\nProcessing domains sequence files")

    # do one more check of pfd files to see if they are all there
    hmm.check_pfd_files(RUN, CLUSTER_BASE_NAMES)

    # BGCs --
    BGCS = big_scape.BGCS() #will contain the BGCs

    if RUN.options.skip_ma:
        print(" Running with skip_ma parameter: Assuming that the domains folder has all the \
            fasta files")
        BGCS.load_from_file(RUN)
    else:
        print(" Adding sequences to corresponding domains file")
        BGCS.load_pfds(RUN, CLUSTER_BASE_NAMES, TRY_RESUME_MULTIPLE_ALIGNMENT)

        BGCS.save_to_file(RUN)

    # Key: BGC. Item: ordered list of simple integers with the number of domains
    # in each gene
    # Instead of `DomainCountGene = defaultdict(list)`, let's try arrays of
    # unsigned ints
    GENE_DOMAIN_COUNT = {}
    # list of gene-numbers that have a hit in the anchor domain list. Zero based
    COREBIOSYNTHETIC_POS = {}
    # list of +/- orientation
    BGC_GENE_ORIENTATION = {}

    # TODO: remove this comment? not sure what it relates to
    # if it's a re-run, the pfd/pfs files were not changed, so the skip_ma flag
    # is activated. We have to open the pfd files to get the gene labels for
    # each domain
    # We now always have to have this data so the alignments are produced
    big_scape.parse_pfd(RUN, CLUSTER_BASE_NAMES, GENE_DOMAIN_COUNT, COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION, BGC_INFO)


    ### Step 5: Create SVG figures
    print(" Creating arrower-like figures for each BGC")

    # read hmm file. We'll need that info anyway for final visualization
    print("  Parsing hmm file for domain information")
    PFAM_INFO = pfam.parse_pfam_a(RUN)
    print("    Done")

    # verify if there are figures already generated

    big_scape.generate_images(RUN, CLUSTER_BASE_NAMES, GEN_BANK_DICT, PFAM_INFO, BGC_INFO)
    print(" Finished creating figures")

    print("\n\n   - - Calculating distance matrix - -")

    # Do multiple alignments if needed
    if not RUN.options.skip_ma:
        hmm.do_multiple_align(RUN, TRY_RESUME_MULTIPLE_ALIGNMENT)

    # If there's something to analyze, load the aligned sequences
    print(" Trying to read domain alignments (*.algn files)")
    ALIGNED_DOMAIN_SEQS = hmm.read_aligned_files(RUN)

    CLUSTER_NAMES = tuple(sorted(CLUSTERS))
    # copy html templates
    dir_util.copy_tree(os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "output"), RUN.directories.output)

    # make a new run folder in the html output & copy the overview_html

    TEMPLATE_PATH = os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "overview_html")
    RUNDATA_NETWORKS_PER_RUN, HTML_SUBS_PER_RUN = big_scape.copy_template_per_cutoff(RUN, TEMPLATE_PATH)

    # create pfams.js
    pfam.create_pfam_js(RUN, PFAM_INFO)

    # Try to make default analysis using all files found inside the input folder
    print("\nGenerating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    print("   Writing the complete Annotations file for the complete set")
    big_scape.write_network_annotation_file(RUN, CLUSTER_NAMES, BGC_INFO)


    # Find index of all MIBiG BGCs if necessary
    if RUN.mibig.use_mibig:
        NAME_TO_IDX = {}
        for clusterIdx, clusterName in enumerate(CLUSTER_NAMES):
            NAME_TO_IDX[clusterName] = clusterIdx

        MIBIG_SET_INDICES = set()
        for bgc in MIBIG_SET:
            MIBIG_SET_INDICES.add(NAME_TO_IDX[bgc])


    # Get the ordered list of domains
    DOMAIN_LIST = big_scape.get_ordered_domain_list(RUN, CLUSTER_BASE_NAMES)

    # Making network files mixing all classes
    if RUN.options.mix:
        big_scape.gen_network_mix_all(RUN, CLUSTER_NAMES, DOMAIN_LIST, BGC_INFO, QUERY_BGC,
                                      GENE_DOMAIN_COUNT, COREBIOSYNTHETIC_POS,
                                      BGC_GENE_ORIENTATION, BGCS, ALIGNED_DOMAIN_SEQS,
                                      MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                      HTML_SUBS_PER_RUN)

    # Making network files separating by BGC class
    if not RUN.options.no_classify:
        big_scape.gen_network_per_class(RUN, CLUSTER_NAMES, DOMAIN_LIST, BGC_INFO, QUERY_BGC,
                                        GENE_DOMAIN_COUNT, COREBIOSYNTHETIC_POS,
                                        BGC_GENE_ORIENTATION, BGCS, ALIGNED_DOMAIN_SEQS,
                                        MIBIG_SET_INDICES, MIBIG_SET, RUNDATA_NETWORKS_PER_RUN,
                                        HTML_SUBS_PER_RUN)

    # fetch genome list for overview.js
    INPUT_CLUSTERS_IDX = []
    big_scape.fetch_genome_list(RUN, INPUT_CLUSTERS_IDX, CLUSTER_NAMES, MIBIG_SET, BGC_INFO, GEN_BANK_DICT)

    # update family data (convert global bgc indexes into input-only indexes)
    big_scape.update_family_data(RUNDATA_NETWORKS_PER_RUN, INPUT_CLUSTERS_IDX, CLUSTER_NAMES, MIBIG_SET)

    # generate overview data
    RUN.end()
    big_scape.generate_results_per_cutoff_value(RUN, RUNDATA_NETWORKS_PER_RUN, HTML_SUBS_PER_RUN)
    # dump bgc info
    pickle.dump(BGC_INFO, open(os.path.join(RUN.directories.cache, 'bgc_info.dict'), 'wb'))

    # done
    RUN.report_runtime()
