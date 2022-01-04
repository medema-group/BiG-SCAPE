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
import sys
import time
from glob import glob
from itertools import combinations
from itertools import product as combinations_product
from collections import defaultdict
from multiprocessing import Pool
import zipfile

from array import array
import json
import shutil
from distutils import dir_util
import networkx as nx
from scipy.optimize.zeros import results_c

# refactored imports
from src import bgctools
from src import big_scape
from src import gbk
from src import hmm
from src import js
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
    OPTIONS = utility.CMD_parser(ROOT_PATH)

    # create new run details
    RUN = big_scape.Run(OPTIONS)

    # run timing starts here
    RUN.start(OPTIONS)

    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    print("\n\n   - - Processing input files - -")

    mibig.extract_mibig(RUN, OPTIONS)

    MIBIG_SET, BGC_INFO, GEN_BANK_DICT = mibig.import_mibig(RUN, OPTIONS)

    print("\nImporting GenBank files")
    gbk.get_gbk_files(RUN.directories.input, RUN, OPTIONS, BGC_INFO, GEN_BANK_DICT, bgctools.BgcData)

    if RUN.directories.has_query_bgc:
        QUERY_BGC = gbk.fileprocessing.import_query_gbk(RUN, OPTIONS, BGC_INFO, GEN_BANK_DICT)
    # clusters and sampleDict contain the necessary structure for all-vs-all and sample analysis
    CLUSTERS = list(GEN_BANK_DICT.keys())

    SAMPLE_DICT = {} # {sampleName:set(bgc1,bgc2,...)}
    GBK_FILES = [] # raw list of gbk file locations
    for (cluster, (path, clusterSample)) in GEN_BANK_DICT.items():
        GBK_FILES.append(path)
        for sample in clusterSample:
            clustersInSample = SAMPLE_DICT.get(sample, set())
            clustersInSample.add(cluster)
            SAMPLE_DICT[sample] = clustersInSample

    print("\nCreating output directories")
    SVG_FOLDER = os.path.join(RUN.directories.output, "SVG")
    utility.create_directory(SVG_FOLDER, "SVG", False)
    NETWORK_FOLDER = os.path.join(RUN.directories.output, "network_files")
    utility.create_directory(NETWORK_FOLDER, "Networks", False)

    print("\nTrying threading on {} cores".format(str(OPTIONS.cores)))

    # BGCs --
    # dictionary of this structure:
    # BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
    #  'specific_domain_name_2'] } }
    # - cluster_name_x: cluster name (can be anything)
    # - general_domain_name_x: PFAM ID, for example 'PF00550'
    # - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in
    # DMS unequivocally. e.g. 'PF00550_start_end', where start and end are genomic positions
    BGCS = {} #will contain the BGCs

    CLASS_NAMES_LEN = len(RUN.distance.bgc_class_names)
    BGC_CLASS_NAME_2_INDEX = dict(zip(RUN.distance.bgc_class_names, range(CLASS_NAMES_LEN)))

    ALIGNED_DOMAIN_SEQS = {} # Key: specific domain sequence label. Item: aligned sequence
    DOMAIN_LIST = {} # Key: BGC. Item: ordered list of domains

    # Key: BGC. Item: ordered list of simple integers with the number of domains
    # in each gene
    # Instead of `DomainCountGene = defaultdict(list)`, let's try arrays of
    # unsigned ints
    GENE_DOMAIN_COUNT = {}
    # list of gene-numbers that have a hit in the anchor domain list. Zero based
    COREBIOSYNTHETIC_POS = {}
    # list of +/- orientation
    BGC_GENE_ORIENTATION = {}

    # to avoid multiple alignment if there's only 1 seq. representing a particular domain
    SEQS_PER_DOMAIN = {}


    ### Step 2: Run hmmscan
    print("\nPredicting domains using hmmscan")

    BASE_NAMES = set(CLUSTERS)

    # All available fasta files (could be more than it should if reusing output folder)
    ALL_FASTA_FILES = set(glob(os.path.join(RUN.directories.bgc_fasta, "*.fasta")))

    # fastaFiles: all the fasta files that should be there
    # (i.e. correspond to the input files)
    FASTA_FILES = set()
    for name in BASE_NAMES:
        FASTA_FILES.add(os.path.join(RUN.directories.bgc_fasta, name+".fasta"))

    # fastaBases: the actual fasta files we have that correspond to the input
    FASTA_BASES = ALL_FASTA_FILES.intersection(FASTA_FILES)

    # Verify that all input files had their fasta sequences extracted
    if len(FASTA_FILES - FASTA_BASES) > 0:
        print("Error! The following files did NOT have their fasta sequences extracted: ")
        UNEXTRACTED_FILES = FASTA_FILES - FASTA_BASES
        for unextracted_file in UNEXTRACTED_FILES:
            print(unextracted_file)
        sys.exit()


    # Make a list of all fasta files that need to be processed
    # (i.e., they don't yet have a corresponding .domtable)
    if OPTIONS.force_hmmscan:
        # process all files, regardless of whether they already existed
        TASK_SET = FASTA_FILES
        print(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
    else:
        # find already processed files
        ALREADY_DONE = set()
        for fasta in FASTA_FILES:
            outputbase = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(RUN.directories.domtable, outputbase + '.domtable')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                # verify domtable content
                with open(outputfile, "r") as domtablefile:
                    for line in domtablefile.readlines():
                        if line.startswith("# Option settings:"):
                            linecols = line.split()
                            if "hmmscan" in linecols and "--domtblout" in linecols:
                                ALREADY_DONE.add(fasta)
                                break

        TASK_SET = FASTA_FILES - ALREADY_DONE
        if len(TASK_SET) == 0:
            print(" All fasta files had already been processed")
        elif len(ALREADY_DONE) > 0:
            if len(TASK_SET) < 20:
                TASKS = [x.split(os.sep)[-1].split(".")[:-1] for x in TASK_SET]
                print(" Warning! The following NEW fasta file(s) will be processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in TASK_SET)))
            else:
                print(" Warning: {} NEW fasta files will be processed".format(len(TASK_SET)))
        else:
            print(" Predicting domains for {} fasta files".format(str(len(FASTA_FILES))))

    POOL = Pool(OPTIONS.cores, maxtasksperchild=1)
    for fastaFile in TASK_SET:
        task_args = (fastaFile, RUN.directories.pfam, RUN.directories.domtable, OPTIONS.verbose)
        POOL.apply_async(hmm.runHmmScan, args=task_args)
    POOL.close()
    POOL.join()
    print(" Finished generating domtable files.")

    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    print("\nParsing hmmscan domtable files")

    # All available domtable files
    ALL_DOMTABLE_FILES = set(glob(os.path.join(RUN.directories.domtable, "*.domtable")))

    # domtableFiles: all domtable files corresponding to the input files
    DOMTABLE_FILES = set()
    for name in BASE_NAMES:
        DOMTABLE_FILES.add(os.path.join(RUN.directories.domtable, name+".domtable"))

    # domtableBases: the actual set of input files with coresponding domtable files
    DOMTABLE_BASES = ALL_DOMTABLE_FILES.intersection(DOMTABLE_FILES)

    # Verify that all input files have a corresponding domtable file
    if len(DOMTABLE_FILES - DOMTABLE_BASES) > 0:
        SKIPPED_FILES = DOMTABLE_FILES - DOMTABLE_BASES
        print("Error! The following files did NOT have their domains predicted: ")
        for skipped_file in SKIPPED_FILES:
            print(skipped_file)
        sys.exit()

    # find already processed files (assuming that if the pfd file exists, the pfs should too)
    ALREADY_DONE = set()
    if not OPTIONS.force_hmmscan:
        for domtable in DOMTABLE_FILES:
            outputbase = ".".join(domtable.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(RUN.directories.pfd, outputbase + '.pfd')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                ALREADY_DONE.add(domtable)
    DOMTABLE_FILES_UNPROCESSED = DOMTABLE_FILES - ALREADY_DONE
    if len(DOMTABLE_FILES_UNPROCESSED) == 0: # Re-run
        print(" All domtable files had already been processed")
    elif len(ALREADY_DONE) > 0: # Incomplete run
        if len(DOMTABLE_FILES_UNPROCESSED) < 20:
            print(" Warning! The following domtable files had not been processed:")
            for unprocessed_domtable_file in DOMTABLE_FILES_UNPROCESSED:
                print(unprocessed_domtable_file.split(os.sep)[-1].split('.')[:-1])
        else:
            print(" Warning: {} domtable files will be processed".format(str(len(DOMTABLE_FILES_UNPROCESSED))))
    else: # First run
        print(" Processing {} domtable files".format(str(len(DOMTABLE_FILES))))

    # If using the multiprocessing version and outputbase doesn't have any
    #  predicted domains, it's not as easy to remove if from the analysis
    #  (probably because parseHmmScan only has a copy of clusters et al?)
    # Using serialized version for now. Probably doesn't have too bad an impact
    #pool = Pool(cores,maxtasksperchild=32)
    for domtableFile in DOMTABLE_FILES - ALREADY_DONE:
        hmm.parseHmmScan(domtableFile, RUN.directories.pfd, RUN.directories.pfs,
                         OPTIONS.domain_overlap_cutoff, OPTIONS.verbose, GEN_BANK_DICT,
                         CLUSTERS, BASE_NAMES, GBK_FILES, SAMPLE_DICT, MIBIG_SET)
        #task_args = (domtableFile,output_folder,options.domain_overlap_cutoff)
        #pool.apply_async(parseHmmScan, args = task_args)
    #pool.close()
    #pool.join()

    # If number of pfd files did not change, no new sequences were added to the
    #  domain fastas and we could try to resume the multiple alignment phase
    # baseNames have been pruned of BGCs with no domains that might've been added temporarily
    TRY_MA_RESUME = False
    if len(BASE_NAMES - set(pfd.split(os.sep)[-1][:-9] for pfd in ALREADY_DONE)) == 0:
        TRY_MA_RESUME = True
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

    # All available pfd files
    ALL_PFD_FILES = set(glob(os.path.join(RUN.directories.pfd, "*.pfd")))

    # pfdFiles: all pfd files corresponding to the input files
    # (some input files could've been removed due to not having predicted domains)
    PFD_FILES = set()
    for name in BASE_NAMES:
        PFD_FILES.add(os.path.join(RUN.directories.pfd, name+".pfd"))

    # pfdBases: the actual set of input files that have pfd files
    PFD_BASES = ALL_PFD_FILES.intersection(PFD_FILES)

    # verify previous step.
    # All BGCs without predicted domains should no longer be in baseNames
    if len(PFD_FILES - PFD_BASES) > 0:
        print("Error! The following files did NOT have their domtable files processed:")
        UNPROCESSED_DOMTABLE_FILES = PFD_FILES - PFD_BASES
        for unprocessed_domtable_file in UNPROCESSED_DOMTABLE_FILES:
            print(unprocessed_domtable_file)
        sys.exit()

    FILTERED_MATRIX = []
    if OPTIONS.skip_ma:
        print(" Running with skip_ma parameter: Assuming that the domains folder has all the fasta \
              files")
        try:
            with open(os.path.join(RUN.directories.cache, "BGCs.dict"), "r") as BGC_file:
                BGCS = pickle.load(BGC_file)
                BGC_file.close()
        except IOError:
            sys.exit("BGCs file not found...")
    else:
        print(" Adding sequences to corresponding domains file")

        for outputbase in BASE_NAMES:
            if OPTIONS.verbose:
                print("   Processing: " + outputbase)

            pfdFile = os.path.join(RUN.directories.pfd, outputbase + ".pfd")
            FILTERED_MATRIX = [[part.strip() for part in l.split('\t')] for l in open(pfdFile)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(RUN.directories.bgc_fasta, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not TRY_MA_RESUME:
                with open(fasta_file, "r") as fasta_file_handle:
                    fasta_dict = utility.fasta_parser(fasta_file_handle) # all fasta info from a BGC
                utility.save_domain_seqs(FILTERED_MATRIX, fasta_dict,
                                         RUN.directories.domains, outputbase)

            BGCS[outputbase] = bgctools.bgc_dict_gen(FILTERED_MATRIX)

            del FILTERED_MATRIX[:]

        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(RUN.directories.cache, "BGCs.dict"), "wb") as BGC_file:
            pickle.dump(BGCS, BGC_file)
            BGC_file.close()

    # if it's a re-run, the pfd/pfs files were not changed, so the skip_ma flag
    # is activated. We have to open the pfd files to get the gene labels for
    # each domain
    # We now always have to have this data so the alignments are produced
    PFD_DICT_DOMAINS = defaultdict(int)
    ORF_KEYS = {}
    for outputbase in BASE_NAMES:
        GENE_DOMAIN_COUNT[outputbase] = array('B')
        COREBIOSYNTHETIC_POS[outputbase] = array('H')
        BGC_GENE_ORIENTATION[outputbase] = array('b')
        pfdFile = os.path.join(RUN.directories.pfd, outputbase + ".pfd")

        #pfd_dict_domains contains the number of domains annotated in the
        # pfd file for each orf tag
        with open(pfdFile, "r") as pfdf:
            for line in pfdf:
                PFD_DICT_DOMAINS[line.strip().split("\t")[-1]] += 1

        # extract the orf number from the tag and use it to traverse the BGC
        for orf in PFD_DICT_DOMAINS.keys():
            orf_num = int(orf.split(":")[0].split("_ORF")[1])
            ORF_KEYS[orf_num] = orf

        orf_num = 0
        for orf_key in sorted(ORF_KEYS.keys()):
            orf = ORF_KEYS[orf_key]
            if orf[-1] == "+":
                BGC_GENE_ORIENTATION[outputbase].append(1)
            else:
                BGC_GENE_ORIENTATION[outputbase].append(-1)

            GENE_DOMAIN_COUNT[outputbase].append(PFD_DICT_DOMAINS[orf])

            if orf in BGC_INFO[outputbase].biosynthetic_genes:
                COREBIOSYNTHETIC_POS[outputbase].append(orf_num)
            orf_num += 1

        PFD_DICT_DOMAINS.clear()
        ORF_KEYS.clear()

        ## TODO: if len(corebiosynthetic_position[outputbase]) == 0
        ## do something with the list of pfam ids. Specifically, mark
        ## (in this case TODO or always?) as biosynthetic genes, the ones that contain
        ## domains from a special list. This list of special domains
        ## comes from predicted domains within the CDSs marked as 'sec_met'
        ## by antismash


    # Get the ordered list of domains
    print(" Reading the ordered list of domains from the pfs files")
    for outputbase in BASE_NAMES:
        pfsfile = os.path.join(RUN.directories.pfs, outputbase + ".pfs")
        if os.path.isfile(pfsfile):
            DOMAIN_LIST[outputbase] = pfam.get_domain_list(pfsfile)
        else:
            sys.exit(" Error: could not open " + outputbase + ".pfs")


    ### Step 5: Create SVG figures
    print(" Creating arrower-like figures for each BGC")

    # read hmm file. We'll need that info anyway for final visualization
    print("  Parsing hmm file for domain information")
    PFAM_INFO = {}
    with open(os.path.join(RUN.directories.pfam, "Pfam-A.hmm"), "r") as pfam:
        PUT_IN_DICT = False
        # assuming that the order of the information never changes
        for line in pfam:
            if line[:4] == "NAME":
                name = line.strip()[6:]
            if line[:3] == "ACC":
                acc = line.strip()[6:].split(".")[0]
            if line[:4] == "DESC":
                desc = line.strip()[6:]
                PUT_IN_DICT = True

            if PUT_IN_DICT:
                PUT_IN_DICT = False
                PFAM_INFO[acc] = (name, desc)
    print("    Done")

    # verify if there are figures already generated

    # All available SVG files
    AVAILABLE_SVG = set()
    for svg in glob(os.path.join(SVG_FOLDER, "*.svg")):
        (root, ext) = os.path.splitext(svg)
        AVAILABLE_SVG.add(root.split(os.sep)[-1])

    # Which files actually need to be generated
    WORKING_SET = BASE_NAMES - AVAILABLE_SVG

    if len(WORKING_SET) > 0:
        COLOR_GENES = {}
        COLOR_DOMAINS = utility.read_color_domains_file()
        PFAM_DOMAIN_CATEGORIES = {}

        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        print("  Reading BGC information and writing SVG")
        for bgc in WORKING_SET:
            with open(GEN_BANK_DICT[bgc][0], "r") as handle:
                utility.SVG(False, os.path.join(SVG_FOLDER, bgc+".svg"), handle, bgc,
                            os.path.join(RUN.directories.pfd, bgc+".pfd"), True, COLOR_GENES,
                            COLOR_DOMAINS, PFAM_DOMAIN_CATEGORIES, PFAM_INFO,
                            BGC_INFO[bgc].records, BGC_INFO[bgc].max_width)

        COLOR_GENES.clear()
        COLOR_DOMAINS.clear()
        PFAM_DOMAIN_CATEGORIES.clear()
    elif len(WORKING_SET) == 0:
        print("  All SVG from the input files seem to be in the SVG folder")

    AVAILABLE_SVG.clear()
    print(" Finished creating figures")


    print("\n\n   - - Calculating distance matrix - -")

    # Do multiple alignments if needed
    if not OPTIONS.skip_ma:
        print("Performing multiple alignment of domain sequences")

        # obtain all fasta files with domain sequences
        DOMAIN_SEQUENCE_LIST = set(glob(os.path.join(RUN.directories.domains, "*.fasta")))

        # compare with .algn set of files. Maybe resuming is possible if
        # no new sequences were added
        if TRY_MA_RESUME:
            TEMP_ALIGNED = set(glob(os.path.join(RUN.directories.domains, "*.algn")))

            if len(TEMP_ALIGNED) > 0:
                print(" Found domain fasta files without corresponding alignments")

                for a in TEMP_ALIGNED:
                    if os.path.getsize(a) > 0:
                        DOMAIN_SEQUENCE_LIST.remove(a[:-5]+".fasta")

            TEMP_ALIGNED.clear()

        # Try to further reduce the set of domain fastas that need alignment
        SEQUENCE_TAG_LIST = set()
        HEADER_LIST = []
        DOMAIN_SEQUENCE_LIST_TEMP = DOMAIN_SEQUENCE_LIST.copy()
        for domain_file in DOMAIN_SEQUENCE_LIST_TEMP:
            domain_name = ".".join(domain_file.split(os.sep)[-1].split(".")[:-1])

            # fill fasta_dict...
            with open(domain_file, "r") as fasta_handle:
                HEADER_LIST = utility.get_fasta_keys(fasta_handle)

            # Get the BGC name from the sequence tag. The form of the tag is:
            # >BGCXXXXXXX_BGCXXXXXXX_ORF25:gid...
            SEQUENCE_TAG_LIST = set(s.split("_ORF")[0] for s in HEADER_LIST)

            # ...to find out how many sequences do we actually have
            if len(SEQUENCE_TAG_LIST) == 1:
                # avoid multiple alignment if the domains all belong to the same BGC
                DOMAIN_SEQUENCE_LIST.remove(domain_file)
                if OPTIONS.verbose:
                    print(" Skipping Multiple Alignment for {} \
                        (appears only in one BGC)".format(domain_name))

        SEQUENCE_TAG_LIST.clear()
        del HEADER_LIST[:]

        DOMAIN_SEQUENCE_LIST_TEMP.clear()

        # Do the multiple alignment
        STOP_FLAG = False
        if len(DOMAIN_SEQUENCE_LIST) > 0:
            print("\n Using hmmalign")
            hmm.launch_hmmalign(OPTIONS.cores, DOMAIN_SEQUENCE_LIST, RUN.directories.pfam,
                                OPTIONS.verbose)

            # verify all tasks were completed by checking existance of alignment files
            for domain_file in DOMAIN_SEQUENCE_LIST:
                if not os.path.isfile(domain_file[:-6]+".algn"):
                    print("   ERROR, {}.algn could not be found \
                        (possible issue with aligner).".format(domain_file[:-6]))
                    STOP_FLAG = True
            if STOP_FLAG:
                sys.exit()
        else:
            print(" No domain fasta files found to align")


    # If there's something to analyze, load the aligned sequences
    print(" Trying to read domain alignments (*.algn files)")
    ALIGNED_FILES_LIST = glob(os.path.join(RUN.directories.domains, "*.algn"))
    if len(ALIGNED_FILES_LIST) == 0:
        sys.exit("No aligned sequences found in the domain folder (run without the --skip_ma parameter or point to the correct output folder)")
    for aligned_file in ALIGNED_FILES_LIST:
        with open(aligned_file, "r") as aligned_file_handle:
            fasta_dict = utility.fasta_parser(aligned_file_handle)
            for header in fasta_dict:
                ALIGNED_DOMAIN_SEQS[header] = fasta_dict[header]

    CLUSTER_NAMES = tuple(sorted(CLUSTERS))

    # we have to find the idx of query_bgc
    if RUN.directories.has_query_bgc:
        try:
            QUERY_BGC_IDX = CLUSTER_NAMES.index(QUERY_BGC)
        except ValueError:
            sys.exit("Error finding the index of Query BGC")

    # create output directory for network files
    NETWORK_FILES_FOLDER = os.path.join(NETWORK_FOLDER, RUN.run_name)
    utility.create_directory(NETWORK_FILES_FOLDER, "Network Files", False)

    # copy html templates
    dir_util.copy_tree(os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "output"), RUN.directories.output)

    # make a new run folder in the html output & copy the overview_html
    NETWORK_HTML_FOLDER = os.path.join(RUN.directories.output, "html_content", "networks", RUN.run_name)
    RUNDATA_NETWORKS_PER_RUN = {}
    HTML_SUBS_PER_RUN = {}
    for cutoff in RUN.cluster.cutoff_list:
        network_html_folder_cutoff = "{}_c{:.2f}".format(NETWORK_HTML_FOLDER, cutoff)
        utility.create_directory(network_html_folder_cutoff, "Network HTML Files", False)
        shutil.copy(os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "overview_html"), os.path.join(network_html_folder_cutoff, "overview.html"))
        RUNDATA_NETWORKS_PER_RUN[network_html_folder_cutoff] = []
        HTML_SUBS_PER_RUN[network_html_folder_cutoff] = []

    # create pfams.js
    PFAMS_JS_FILE = os.path.join(RUN.directories.output, "html_content", "js", "pfams.js")
    if not os.path.isfile(PFAMS_JS_FILE):
        with open(PFAMS_JS_FILE, "w") as pfams_js:
            PFAM_JSON = {}
            PFAM_COLORS = pfam.generatePfamColorsMatrix(os.path.join(os.path.dirname(os.path.realpath(__file__)), "domains_color_file.tsv"))
            for pfam_code in PFAM_INFO:
                pfam_obj = {}
                if pfam_code in PFAM_COLORS:
                    pfam_obj["col"] = PFAM_COLORS[pfam_code]
                else:
                    pfam_obj["col"] = "255,255,255"
                pfam_obj["desc"] = PFAM_INFO[pfam_code][1]
                PFAM_JSON[pfam_code] = pfam_obj
            pfams_js.write("var pfams={};\n".format(json.dumps(PFAM_JSON, indent=4, separators=(',', ':'), sort_keys=True)))

    # Try to make default analysis using all files found inside the input folder
    print("\nGenerating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    print("   Writing the complete Annotations file for the complete set")
    NETWORK_ANNOTATION_PATH = os.path.join(NETWORK_FILES_FOLDER, "Network_Annotations_Full.tsv")
    with open(NETWORK_ANNOTATION_PATH, "w") as network_annotation_file:
        network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
        for bgc in CLUSTER_NAMES:
            product = BGC_INFO[bgc].product
            network_annotation_file.write("\t".join([bgc, BGC_INFO[bgc].accession_id, BGC_INFO[bgc].description, product, bgctools.sort_bgc(product), BGC_INFO[bgc].organism, BGC_INFO[bgc].taxonomy]) + "\n")


    # Find index of all MIBiG BGCs if necessary
    if RUN.mibig.use_mibig:
        NAME_TO_IDX = {}
        for clusterIdx, clusterName in enumerate(CLUSTER_NAMES):
            NAME_TO_IDX[clusterName] = clusterIdx

        MIBIG_SET_INDICES = set()
        for bgc in MIBIG_SET:
            MIBIG_SET_INDICES.add(NAME_TO_IDX[bgc])

    # Making network files mixing all classes
    if OPTIONS.mix:
        print("\n Mixing all BGC classes")

        # only choose from valid classes
        MIX_SET = []

        # create working set with indices of valid clusters
        for clusterIdx, clusterName in enumerate(CLUSTER_NAMES):
            if RUN.has_includelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DOMAIN_LIST[clusterName]})

                if len(RUN.domain_includelist & bgc_domain_set) == 0:
                    continue

            product = BGC_INFO[clusterName].product
            predicted_class = bgctools.sort_bgc(product)

            if predicted_class.lower() in RUN.valid_classes:
                MIX_SET.append(clusterIdx)

        print("\n  {} ({} BGCs)".format("Mix", str(len(MIX_SET))))

        # create output directory
        utility.create_directory(os.path.join(NETWORK_FILES_FOLDER, "mix"), "  Mix", False)

        print("  Calculating all pairwise distances")
        if RUN.directories.has_query_bgc:
            PAIRS = set([tuple(sorted(combo)) for combo in combinations_product([QUERY_BGC_IDX], MIX_SET)])
        else:
            # convert into a set of ordered tuples
            PAIRS = set([tuple(sorted(combo)) for combo in combinations(MIX_SET, 2)])

        CLUSTER_PAIRS = [(x, y, -1) for (x, y) in PAIRS]
        PAIRS.clear()
        NETWORK_MATRIX_MIX = big_scape.generate_network(CLUSTER_PAIRS, OPTIONS.cores,
                                                        CLUSTER_NAMES,
                                                        RUN.distance.bgc_class_names,
                                                        DOMAIN_LIST, RUN.directories.output,
                                                        GENE_DOMAIN_COUNT,
                                                        COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION,
                                                        RUN.distance.bgc_class_weight,
                                                        RUN.network.anchor_domains, BGCS,
                                                        OPTIONS.mode, BGC_INFO,
                                                        ALIGNED_DOMAIN_SEQS, OPTIONS.verbose,
                                                        RUN.directories.domains)

        del CLUSTER_PAIRS[:]

        # add the rest of the edges in the "Query network"
        if RUN.directories.has_query_bgc:
            NEW_SET = []

            # rows from the distance matrix that will be pruned
            DEL_LIST = []

            for idx, row in enumerate(NETWORK_MATRIX_MIX):
                a, b, distance = int(row[0]), int(row[1]), row[2]

                if a == b:
                    continue

                if distance <= RUN.cluster.max_cutoff:
                    if a == QUERY_BGC_IDX:
                        NEW_SET.append(b)
                    else:
                        NEW_SET.append(a)
                else:
                    DEL_LIST.append(idx)

            for idx in sorted(DEL_LIST, reverse=True):
                del NETWORK_MATRIX_MIX[idx]
            del DEL_LIST[:]

            PAIRS = set([tuple(sorted(combo)) for combo in combinations(NEW_SET, 2)])
            CLUSTER_PAIRS = [(x, y, -1) for (x, y) in PAIRS]
            PAIRS.clear()
            NETWORK_MATRIX_NEW_SET = big_scape.generate_network(CLUSTER_PAIRS, OPTIONS.cores, CLUSTER_NAMES, RUN.distance.bgc_class_names, DOMAIN_LIST, RUN.directories.output, GENE_DOMAIN_COUNT,
            COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION, RUN.distance.bgc_class_weight, RUN.network.anchor_domains, BGCS, OPTIONS.mode, BGC_INFO,
            ALIGNED_DOMAIN_SEQS, OPTIONS.verbose, RUN.directories.domains)
            del CLUSTER_PAIRS[:]

            # Update the network matrix (QBGC-vs-all) with the distances of
            # QBGC's GCF
            NETWORK_MATRIX_MIX.extend(NETWORK_MATRIX_NEW_SET)

            # Update actual list of BGCs that we'll use
            MIX_SET = NEW_SET
            MIX_SET.extend([QUERY_BGC_IDX])
            MIX_SET.sort() # clusterJsonBatch expects ordered indices

            # Create an additional file with the list of all clusters in the class + other info
            # This version of the file only has information on the BGCs connected to Query BGC
            print("   Writing annotation file")
            NETWORK_ANNOTATION_PATH = os.path.join(NETWORK_FILES_FOLDER, "mix", "Network_Annotations_mix_QueryBGC.tsv")
            with open(NETWORK_ANNOTATION_PATH, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in MIX_SET:
                    bgc = CLUSTER_NAMES[idx]
                    product = BGC_INFO[bgc].product
                    network_annotation_file.write("\t".join([bgc,
                        BGC_INFO[bgc].accession_id, BGC_INFO[bgc].description,
                        product, bgctools.sort_bgc(product), BGC_INFO[bgc].organism,
                        BGC_INFO[bgc].taxonomy]) + "\n")
        elif RUN.mibig.use_mibig:
            N = nx.Graph()
            N.add_nodes_from(MIX_SET)
            MIBIG_SET_DEL = []
            NETWORK_MATRIX_SET_DEL = []

            for idx, row in enumerate(NETWORK_MATRIX_MIX):
                a, b, distance = int(row[0]), int(row[1]), row[2]
                if distance <= RUN.cluster.max_cutoff:
                    N.add_edge(a, b, index=idx)

            for component in nx.connected_components(N): # note: 'component' is a set
                numBGCs_subgraph = len(component)

                # catch if the subnetwork is comprised only of MIBiG BGCs
                if len(component & MIBIG_SET_INDICES) == numBGCs_subgraph:
                    for bgc in component:
                        MIBIG_SET_DEL.append(bgc)

            # Get all edges between bgcs marked for deletion
            for (a, b, idx) in N.subgraph(MIBIG_SET_DEL).edges.data('index'):
                NETWORK_MATRIX_SET_DEL.append(idx)

            # delete all edges between marked bgcs
            for row_idx in sorted(NETWORK_MATRIX_SET_DEL, reverse=True):
                del NETWORK_MATRIX_MIX[row_idx]
            del NETWORK_MATRIX_SET_DEL[:]

            print("   Removing {} non-relevant MIBiG BGCs".format(len(MIBIG_SET_DEL)))
            MIX_SET_IDX = 0
            BGC_TO_MIX_SET_IDX = {}
            for idx, bgc in enumerate(MIX_SET):
                BGC_TO_MIX_SET_IDX[bgc] = idx

            for bgc_idx in sorted(MIBIG_SET_DEL, reverse=True):
                del MIX_SET[BGC_TO_MIX_SET_IDX[bgc_idx]]
            del MIBIG_SET_DEL[:]


        print("  Writing output files")
        PATH_BASE = os.path.join(NETWORK_FILES_FOLDER, "mix")
        FILE_NAMES = []
        for cutoff in RUN.cluster.cutoff_list:
            FILE_NAMES.append(os.path.join(PATH_BASE, "mix_c{:.2f}.network".format(cutoff)))
        CUTOFFS_FILENAMES = list(zip(RUN.cluster.cutoff_list, FILE_NAMES))
        del FILE_NAMES[:]
        big_scape.write_network_matrix(NETWORK_MATRIX_MIX, CUTOFFS_FILENAMES, OPTIONS.include_singletons, CLUSTER_NAMES, BGC_INFO)

        print("  Calling Gene Cluster Families")
        REDUCED_NETWORK = []
        POS_ALIGNMENTS = {}
        for row in NETWORK_MATRIX_MIX:
            REDUCED_NETWORK.append([int(row[0]), int(row[1]), row[2]])
            reverse = False
            if row[-1] == 1.0:
                reverse = True
            pa = POS_ALIGNMENTS.setdefault(int(row[0]), {})
            # lcsStartA, lcsStartB, seedLength, reverse={True,False}
            pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
        del NETWORK_MATRIX_MIX[:]
        FAMILY_DATA = big_scape.clusterJsonBatch(MIX_SET, PATH_BASE, "mix", REDUCED_NETWORK, POS_ALIGNMENTS,
            CLUSTER_NAMES, BGC_INFO, MIBIG_SET, RUN.directories.pfd, RUN.directories.bgc_fasta,
            DOMAIN_LIST, BGCS, ALIGNED_DOMAIN_SEQS, GENE_DOMAIN_COUNT, BGC_GENE_ORIENTATION,
            cutoffs=RUN.cluster.cutoff_list, clusterClans=OPTIONS.clans,
            clanCutoff=OPTIONS.clan_cutoff, htmlFolder=NETWORK_HTML_FOLDER)
        for network_html_folder_cutoff in FAMILY_DATA:
            RUNDATA_NETWORKS_PER_RUN[network_html_folder_cutoff].append(FAMILY_DATA[network_html_folder_cutoff])
            HTML_SUBS_PER_RUN[network_html_folder_cutoff].append({"name" : "mix", "css" : "Others", "label" : "Mixed"})
        del MIX_SET[:]
        del REDUCED_NETWORK[:]


    # Making network files separating by BGC class
    if not OPTIONS.no_classify:
        print("\n Working for each BGC class")

        # TODO: remove?
        # reinitialize run.distance.bgc_classes to make sure the bgc lists are empty
        RUN.distance.bgc_classes = defaultdict(list)

        # Preparing gene cluster classes
        print("  Sorting the input BGCs\n")

        # create and sort working set for each class
        for clusterIdx, clusterName in enumerate(CLUSTER_NAMES):
            if RUN.has_includelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DOMAIN_LIST[clusterName]})

                if len(RUN.domain_includelist & bgc_domain_set) == 0:
                    continue

            product = BGC_INFO[clusterName].product
            predicted_class = bgctools.sort_bgc(product)

            if predicted_class.lower() in RUN.valid_classes:
                RUN.distance.bgc_classes[predicted_class].append(clusterIdx)

            # possibly add hybrids to 'pure' classes
            if OPTIONS.hybrids:
                if predicted_class == "PKS-NRP_Hybrids":
                    if "nrps" in RUN.valid_classes:
                        RUN.distance.bgc_classes["NRPS"].append(clusterIdx)
                    if "t1pks" in product and "pksi" in RUN.valid_classes:
                        RUN.distance.bgc_classes["PKSI"].append(clusterIdx)
                    if "t1pks" not in product and "pksother" in RUN.valid_classes:
                        RUN.distance.bgc_classes["PKSother"].append(clusterIdx)

                if predicted_class == "Others" and "." in product:
                    subclasses = set()
                    for subproduct in product.split("."):
                        subclass = bgctools.sort_bgc(subproduct)
                        if subclass.lower() in RUN.valid_classes:
                            subclasses.add(subclass)

                    # Prevent mixed BGCs with sub-Others annotations to get
                    # added twice (e.g. indole-cf_fatty_acid has already gone
                    # to Others at this point)
                    if "Others" in subclasses:
                        subclasses.remove("Others")


                    for subclass in subclasses:
                        RUN.distance.bgc_classes[subclass].append(clusterIdx)
                    subclasses.clear()

        # only make folders for the run.distance.bgc_classes that are found
        for bgc_class in RUN.distance.bgc_classes:
            if RUN.directories.has_query_bgc:
                # not interested in this class if our Query BGC is not here...
                if QUERY_BGC_IDX not in RUN.distance.bgc_classes[bgc_class]:
                    continue

            print("\n  {} ({} BGCs)".format(bgc_class, str(len(RUN.distance.bgc_classes[bgc_class]))))
            if RUN.mibig.use_mibig:
                if len(set(RUN.distance.bgc_classes[bgc_class]) & MIBIG_SET_INDICES) == len(RUN.distance.bgc_classes[bgc_class]):
                    print(" - All clusters in this class are MIBiG clusters -")
                    print("  If you'd like to analyze MIBiG clusters, turn off the --mibig option")
                    print("  and point --inputdir to the Annotated_MIBiG_reference folder")
                    continue

            # create output directory
            utility.create_directory(os.path.join(NETWORK_FILES_FOLDER, bgc_class), "  All - " + bgc_class, False)

            # Create an additional file with the final list of all clusters in the class
            print("   Writing annotation files")
            NETWORK_ANNOTATION_PATH = os.path.join(NETWORK_FILES_FOLDER, bgc_class, "Network_Annotations_" + bgc_class + ".tsv")
            with open(NETWORK_ANNOTATION_PATH, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in RUN.distance.bgc_classes[bgc_class]:
                    bgc = CLUSTER_NAMES[idx]
                    product = BGC_INFO[bgc].product
                    network_annotation_file.write("\t".join([bgc, BGC_INFO[bgc].accession_id, BGC_INFO[bgc].description, product, bgctools.sort_bgc(product), BGC_INFO[bgc].organism, BGC_INFO[bgc].taxonomy]) + "\n")

            print("   Calculating all pairwise distances")
            if RUN.directories.has_query_bgc:
                PAIRS = set([tuple(sorted(combo)) for combo in combinations_product([QUERY_BGC_IDX], RUN.distance.bgc_classes[bgc_class])])
            else:
                PAIRS = set([tuple(sorted(combo)) for combo in combinations(RUN.distance.bgc_classes[bgc_class], 2)])

            CLUSTER_PAIRS = [(x, y, BGC_CLASS_NAME_2_INDEX[bgc_class]) for (x, y) in PAIRS]
            PAIRS.clear()
            network_matrix = big_scape.generate_network(CLUSTER_PAIRS, OPTIONS.cores, CLUSTER_NAMES, RUN.distance.bgc_class_names, DOMAIN_LIST, RUN.directories.output, GENE_DOMAIN_COUNT,
            COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION, RUN.distance.bgc_class_weight, RUN.network.anchor_domains, BGCS, OPTIONS.mode, BGC_INFO,
            ALIGNED_DOMAIN_SEQS, OPTIONS.verbose, RUN.directories.domains)
            #pickle.dump(network_matrix,open("others.ntwrk",'wb'))
            del CLUSTER_PAIRS[:]
            #network_matrix = pickle.load(open("others.ntwrk", "rb"))

            # add the rest of the edges in the "Query network"
            if RUN.directories.has_query_bgc:
                NEW_SET = []

                # rows from the distance matrix that will be pruned
                DEL_LIST = []

                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]

                    # avoid QBGC-QBGC
                    if a == b:
                        continue

                    if distance <= RUN.cluster.max_cutoff:
                        if a == QUERY_BGC_IDX:
                            NEW_SET.append(b)
                        else:
                            NEW_SET.append(a)
                    else:
                        DEL_LIST.append(idx)

                for idx in sorted(DEL_LIST, reverse=True):
                    del network_matrix[idx]
                del DEL_LIST[:]

                PAIRS = set([tuple(sorted(combo)) for combo in combinations(NEW_SET, 2)])
                CLUSTER_PAIRS = [(x, y, BGC_CLASS_NAME_2_INDEX[bgc_class]) for (x, y) in PAIRS]
                PAIRS.clear()
                NETWORK_MATRIX_NEW_SET = big_scape.generate_network(CLUSTER_PAIRS, OPTIONS.cores, CLUSTER_NAMES, RUN.distance.bgc_class_names, DOMAIN_LIST, RUN.directories.output, GENE_DOMAIN_COUNT,
                COREBIOSYNTHETIC_POS, BGC_GENE_ORIENTATION, RUN.distance.bgc_class_weight, RUN.network.anchor_domains, BGCS, OPTIONS.mode, BGC_INFO,
                ALIGNED_DOMAIN_SEQS, OPTIONS.verbose, RUN.directories.domains)
                del CLUSTER_PAIRS[:]

                # Update the network matrix (QBGC-vs-all) with the distances of
                # QBGC's GCF
                network_matrix.extend(NETWORK_MATRIX_NEW_SET)

                # Update actual list of BGCs that we'll use
                RUN.distance.bgc_classes[bgc_class] = NEW_SET
                RUN.distance.bgc_classes[bgc_class].extend([QUERY_BGC_IDX])
                RUN.distance.bgc_classes[bgc_class].sort()

                # Create an additional file with the list of all clusters in the class + other info
                # This version of the file only has information on the BGCs connected to Query BGC
                print("   Writing annotation file (Query BGC)")
                NETWORK_ANNOTATION_PATH = os.path.join(NETWORK_FILES_FOLDER, bgc_class, "Network_Annotations_" + bgc_class + "_QueryBGC.tsv")
                with open(NETWORK_ANNOTATION_PATH, "w") as network_annotation_file:
                    network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                    for idx in RUN.distance.bgc_classes[bgc_class]:
                        bgc = CLUSTER_NAMES[idx]
                        product = BGC_INFO[bgc].product
                        network_annotation_file.write("\t".join([bgc, BGC_INFO[bgc].accession_id, BGC_INFO[bgc].description, product, bgctools.sort_bgc(product), BGC_INFO[bgc].organism, BGC_INFO[bgc].taxonomy]) + "\n")
            elif RUN.mibig.use_mibig:
                N = nx.Graph()
                N.add_nodes_from(RUN.distance.bgc_classes[bgc_class])
                MIBIG_SET_DEL = []
                NETWORK_MATRIX_SET_DEL = []

                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]
                    if distance <= RUN.cluster.max_cutoff:
                        N.add_edge(a, b, index=idx)

                for component in nx.connected_components(N): # note: 'component' is a set
                    numBGCs_subgraph = len(component)

                    # catch if the subnetwork is comprised only of MIBiG BGCs
                    if len(component & MIBIG_SET_INDICES) == numBGCs_subgraph:
                        for bgc in component:
                            MIBIG_SET_DEL.append(bgc)

                # Get all edges between bgcs marked for deletion
                for (a, b, idx) in N.subgraph(MIBIG_SET_DEL).edges.data('index'):
                    NETWORK_MATRIX_SET_DEL.append(idx)

                # delete all edges between marked bgcs
                for row_idx in sorted(NETWORK_MATRIX_SET_DEL, reverse=True):
                    del network_matrix[row_idx]
                del NETWORK_MATRIX_SET_DEL[:]

                print("   Removing {} non-relevant MIBiG BGCs".format(len(MIBIG_SET_DEL)))
                bgc_to_class_idx = {}
                for idx, bgc in enumerate(RUN.distance.bgc_classes[bgc_class]):
                    bgc_to_class_idx[bgc] = idx
                for bgc_idx in sorted(MIBIG_SET_DEL, reverse=True):
                    del RUN.distance.bgc_classes[bgc_class][bgc_to_class_idx[bgc_idx]]
                del MIBIG_SET_DEL[:]



            if len(RUN.distance.bgc_classes[bgc_class]) < 2:
                continue

            print("   Writing output files")
            PATH_BASE = os.path.join(NETWORK_FILES_FOLDER, bgc_class)
            FILE_NAMES = []
            for cutoff in RUN.cluster.cutoff_list:
                FILE_NAMES.append(os.path.join(PATH_BASE, "{}_c{:.2f}.network".format(bgc_class, cutoff)))
            CUTOFFS_FILENAMES = list(zip(RUN.cluster.cutoff_list, FILE_NAMES))
            del FILE_NAMES[:]
            big_scape.write_network_matrix(network_matrix, CUTOFFS_FILENAMES, OPTIONS.include_singletons, CLUSTER_NAMES, BGC_INFO)

            print("  Calling Gene Cluster Families")
            REDUCED_NETWORK = []
            POS_ALIGNMENTS = {}
            for row in network_matrix:
                REDUCED_NETWORK.append([int(row[0]), int(row[1]), row[2]])
                reverse = False
                if row[-1] == 1.0:
                    reverse = True
                pa = POS_ALIGNMENTS.setdefault(int(row[0]), {})
                # lcsStartA, lcsStartB, seedLength, reverse={True,False}
                pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
            del network_matrix[:]

            FAMILY_DATA = big_scape.clusterJsonBatch(RUN.distance.bgc_classes[bgc_class], PATH_BASE, bgc_class,
                REDUCED_NETWORK, POS_ALIGNMENTS, CLUSTER_NAMES, BGC_INFO,
                MIBIG_SET, RUN.directories.pfd, RUN.directories.bgc_fasta, DOMAIN_LIST,
                BGCS, ALIGNED_DOMAIN_SEQS, GENE_DOMAIN_COUNT, BGC_GENE_ORIENTATION,
                cutoffs=RUN.cluster.cutoff_list, clusterClans=OPTIONS.clans, clanCutoff=OPTIONS.clan_cutoff,
                htmlFolder=NETWORK_HTML_FOLDER)
            for network_html_folder_cutoff in FAMILY_DATA:
                RUNDATA_NETWORKS_PER_RUN[network_html_folder_cutoff].append(FAMILY_DATA[network_html_folder_cutoff])
                if len(FAMILY_DATA[network_html_folder_cutoff]["families"]) > 0:
                    HTML_SUBS_PER_RUN[network_html_folder_cutoff].append({"name" : bgc_class, "css" : bgc_class, "label" : bgc_class})
            del RUN.distance.bgc_classes[bgc_class][:]
            del REDUCED_NETWORK[:]

    # fetch genome list for overview.js
    GENOMES = []
    CLASSES = []
    CLUSTER_NAMES_TO_GENOMES = {}
    CLUSTER_NAMES_TO_CLASSES = {}
    INPUT_CLUSTERS_IDX = [] # contain only indexes (from clusterNames) of input BGCs (non-mibig)
    for idx, bgc in enumerate(CLUSTER_NAMES):
        if bgc in MIBIG_SET:
            continue
        INPUT_CLUSTERS_IDX.append(idx)
        # get class info
        product = BGC_INFO[bgc].product
        predicted_class = bgctools.sort_bgc(product)
        if predicted_class not in CLASSES:
            CLUSTER_NAMES_TO_CLASSES[bgc] = len(CLASSES)
            CLASSES.append(predicted_class)
        else:
            CLUSTER_NAMES_TO_CLASSES[bgc] = CLASSES.index(predicted_class)
        # get identifier info
        identifier = ""
        if len(BGC_INFO[bgc].organism) > 1:
            identifier = BGC_INFO[bgc].organism
        else: # use original genome file name (i.e. exclude "..clusterXXX from antiSMASH run")
            file_name_base = os.path.splitext(os.path.basename(GEN_BANK_DICT[bgc][0]))[0]
            identifier = file_name_base.rsplit(".cluster", 1)[0].rsplit(".region", 1)[0]
        if len(identifier) < 1:
            identifier = "Unknown Genome {}".format(len(GENOMES))
        if identifier not in GENOMES:
            CLUSTER_NAMES_TO_GENOMES[bgc] = len(GENOMES)
            GENOMES.append(identifier)
        else:
            CLUSTER_NAMES_TO_GENOMES[bgc] = GENOMES.index(identifier)
    RUN.run_data["input"]["accession"] = [{"id": "genome_{}".format(i), "label": acc} for i, acc in enumerate(GENOMES)]
    RUN.run_data["input"]["accession_newick"] = [] # todo ...
    RUN.run_data["input"]["classes"] = [{"label": cl} for cl in CLASSES] # todo : colors
    RUN.run_data["input"]["bgc"] = [{"id": CLUSTER_NAMES[idx], "acc": CLUSTER_NAMES_TO_GENOMES[CLUSTER_NAMES[idx]], "class": CLUSTER_NAMES_TO_CLASSES[CLUSTER_NAMES[idx]]} for idx in INPUT_CLUSTERS_IDX]


    # update family data (convert global bgc indexes into input-only indexes)
    for network_key in RUNDATA_NETWORKS_PER_RUN:
        for network in RUNDATA_NETWORKS_PER_RUN[network_key]:
            for family in network["families"]:
                new_members = []
                mibig = []
                for bgcIdx in family["members"]:
                    if bgcIdx in INPUT_CLUSTERS_IDX:
                        new_members.append(INPUT_CLUSTERS_IDX.index(bgcIdx))
                    else: # is a mibig bgc
                        clusterName = CLUSTER_NAMES[bgcIdx]
                        if clusterName in MIBIG_SET:
                            mibig.append(clusterName)
                family["mibig"] = mibig
                family["members"] = new_members


    # generate overview data
    END_TIME = time.time()
    DURATION = END_TIME - RUN.start_time
    RUN.run_data["end_time"] = time.strftime("%d/%m/%Y %H:%M:%S", time.localtime(END_TIME))
    RUN.run_data["duration"] = "{}h{}m{}s".format((DURATION // 3600), ((DURATION % 3600) // 60), ((DURATION % 3600) % 60))

    for cutoff in RUN.cluster.cutoff_list:
        # update overview.html
        html_folder_for_this_cutoff = "{}_c{:.2f}".format(NETWORK_HTML_FOLDER, cutoff)
        run_data_for_this_cutoff = RUN.run_data.copy()
        run_data_for_this_cutoff["networks"] = RUNDATA_NETWORKS_PER_RUN[html_folder_for_this_cutoff]
        with open(os.path.join(html_folder_for_this_cutoff, "run_data.js"), "w") as run_data_js:
            run_data_js.write("var run_data={};\n".format(json.dumps(run_data_for_this_cutoff, indent=4, separators=(',', ':'), sort_keys=True)))
            run_data_js.write("dataLoaded();\n")
        # update bgc_results.js
        RUN_STRING = "{}_c{:.2f}".format(RUN.run_name, cutoff)
        RESULTS_PATH = os.path.join(RUN.directories.output, "html_content", "js",
                                    "bigscape_results.js")
        js.add_to_bigscape_results_js(RUN_STRING, HTML_SUBS_PER_RUN[html_folder_for_this_cutoff],
                                      RESULTS_PATH)

    pickle.dump(BGC_INFO, open(os.path.join(RUN.directories.cache, 'bgc_info.dict'), 'wb'))
    RUNTIME = time.time()-RUN.start_time
    RUNTIME_STRING = "\n\n\tMain function took {:.3f} s".format(RUNTIME)
    with open(os.path.join(RUN.directories.log, "runtimes.txt"), 'a') as timings_file:
        timings_file.write(RUNTIME_STRING + "\n")
    print(RUNTIME_STRING)
