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

from sys import version_info
if version_info[0]==2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0]==3:
    import pickle # for storing and retrieving dictionaries

import os
import sys
import time
from glob import glob
from itertools import combinations
from itertools import product as combinations_product
from collections import defaultdict
from multiprocessing import Pool
import zipfile

from src.utility.ArrowerSVG import *

from array import array
import json
import shutil
from distutils import dir_util
import networkx as nx



# refactored imports
from src.fileprocessing.gbk import get_gbk_files

from src.utility.cmd_parser import CMD_parser
from src.utility.misc import get_anchor_domains
from src.utility.io import create_directory, write_parameters
from src.utility.fasta import fasta_parser, save_domain_seqs, get_fasta_keys

from src.bgc.misc import sort_bgc, BGC_dic_gen, bgc_data

from src.hmm.hmmscan import runHmmScan, parseHmmScan
from src.hmm.hmmalign import launch_hmmalign

from src.pfam.misc import get_domain_list, generatePfamColorsMatrix

from src.big_scape.network import generate_network, write_network_matrix
from src.big_scape.clustering import clusterJsonBatch

from src.js.misc import add_to_bigscape_results_js




if __name__=="__main__":
    options = CMD_parser()
    
    if options.outputdir == "":
        print("please provide a name for an output folder using parameter -o or --outputdir")
        sys.exit(0)
    
    if os.path.isfile(options.anchorfile):
        anchor_domains = get_anchor_domains(options.anchorfile)
    else:
        print("File with list of anchor domains not found")
        anchor_domains = set()
    
    
    include_singletons = options.include_singletons
    
    cores = int(options.cores)
    
    options_mix = options.mix
    options_classify = not options.no_classify
    
    force_hmmscan = options.force_hmmscan
    
    mode = options.mode
    
    cutoff_list = options.cutoffs
    for c in cutoff_list:
        if c <= 0.0 or c > 1.0:
            sys.exit(" Invalid cutoff value {}".format(str(c)))
    max_cutoff = max(cutoff_list)
            
    # if we want to classify by clans make sure that the clanCutoff is included
    # in the cutoffs to do AP in clusterJsonBatch
    if options.clans:
        fc, cc = options.clan_cutoff
        if cc < fc:
            sys.exit("Error: first value in the clan_cutoff parameter should be smaller than the second")
        if fc not in cutoff_list:
            if fc <= 0.0 or fc > 1.0:
                sys.exit("Error: invalid cutoff value for GCF calling")
            else:
                cutoff_list.append(fc)
                cutoff_list.sort()
            
        if cc <= 0.0 or cc > 1.0:
            sys.exit("Error: invalid cutoff value for GCC calling")
        
    output_folder = str(options.outputdir)
    
    pfam_dir = str(options.pfam_dir)
    h3f = os.path.join(pfam_dir, "Pfam-A.hmm.h3f")
    h3i = os.path.join(pfam_dir, "Pfam-A.hmm.h3i")
    h3m = os.path.join(pfam_dir, "Pfam-A.hmm.h3m")
    h3p = os.path.join(pfam_dir, "Pfam-A.hmm.h3p")
    if not (os.path.isfile(h3f) and os.path.isfile(h3i) and os.path.isfile(h3m) and os.path.isfile(h3p)):
        print("One or more of the necessary Pfam files (.h3f, .h3i, .h3m, .h3p) were not found")
        if os.path.isfile(os.path.join(pfam_dir, "Pfam-A.hmm")):
            print("Please use hmmpress with Pfam-A.hmm")
        else:
            print("Please download the latest Pfam-A.hmm file from http://pfam.xfam.org/")
            print("Then use hmmpress on it, and use the --pfam_dir parameter to point to the location of the files")
        sys.exit()

    has_query_bgc = False
    if options.query_bgc:
        has_query_bgc = True
        if not os.path.isfile(options.query_bgc):
            sys.exit("Error: Query BGC not found")
            
    verbose = options.verbose
    
    selected_mibig = 0
    if options.mibig21: selected_mibig += 1
    if options.mibig14: selected_mibig += 1
    if options.mibig13: selected_mibig += 1
    if selected_mibig > 1:
        sys.exit("Error: choose only one MIBiG version")
    use_relevant_mibig = False
    if selected_mibig == 1:
        use_relevant_mibig = True
    
    run_mode_string = ""
    networks_folder_all = "networks_all"
    if options.hybrids:
        networks_folder_all += "_hybrids"
        run_mode_string += "_hybrids"
    if mode == "auto":
        networks_folder_all += "_auto"
        run_mode_string += "_auto"
    elif mode == "glocal":
        networks_folder_all += "_glocal"
        run_mode_string += "_glocal"
    else:
        run_mode_string += "_global"
    
    time1 = time.time()

    start_time = time.localtime()
    run_name = "{}{}".format(time.strftime("%Y-%m-%d_%H-%M-%S", start_time), run_mode_string)
    if options.label:
        run_name = run_name + "_" + options.label
    run_data = {}
    run_data["start_time"] = time.strftime("%d/%m/%Y %H:%M:%S", start_time)
    run_data["parameters"] = " ".join(sys.argv[1:])
    run_data["input"] = {}
    
    
    # Get domain_includelist
    has_includelist = False
    if options.domain_includelist:
        bigscape_path = os.path.dirname(os.path.realpath(__file__))
        if os.path.isfile(os.path.join(bigscape_path,"domain_includelist.txt")):
            domain_includelist = set()
            for line in open(os.path.join(bigscape_path,"domain_includelist.txt"), "r"):
                if line[0] == "#":
                    continue
                domain_includelist.add(line.split("\t")[0].strip())
            if len(domain_includelist) == 0:
                print("Warning: --domain_includelist used, but no domains found in the file")
            else:
                has_includelist = True
        else:
            sys.exit("Error: domain_includelist.txt file not found")
    
    
    ### Step 1: Get all the input files. Write extract sequence and write fasta if necessary
    print("\n\n   - - Processing input files - -")
    
    create_directory(output_folder, "Output", False)

    # logs
    log_folder = os.path.join(output_folder, "logs")
    create_directory(log_folder, "Logs", False)
    write_parameters(log_folder, sys.argv)

    # cached stuff
    cache_folder = os.path.join(output_folder, "cache")
    bgc_fasta_folder = os.path.join(cache_folder, "fasta")
    domtable_folder = os.path.join(cache_folder, "domtable")
    pfs_folder = os.path.join(cache_folder, "pfs")
    pfd_folder = os.path.join(cache_folder, "pfd")
    domains_folder = os.path.join(cache_folder, "domains")
    create_directory(cache_folder, "Cache", False)
    create_directory(bgc_fasta_folder, "BGC fastas", False)
    create_directory(domtable_folder, "Domtable", False)
    create_directory(domains_folder, "Domains", False)
    create_directory(pfs_folder, "pfs", False)
    create_directory(pfd_folder, "pfd", False)
    
    
    # Weights in the format J, DSS, AI, anchorboost
    # Generated with optimization results 2016-12-05. 
    # Used the basic list of 4 anchor domains.
    bgc_class_weight = {}
    bgc_class_weight["PKSI"] = (0.22, 0.76, 0.02, 1.0)
    bgc_class_weight["PKSother"] = (0.0, 0.32, 0.68, 4.0)
    bgc_class_weight["NRPS"] = (0.0, 1.0, 0.0, 4.0)
    bgc_class_weight["RiPPs"] = (0.28, 0.71, 0.01, 1.0)
    bgc_class_weight["Saccharides"] = (0.0, 0.0, 1.0, 1.0)
    bgc_class_weight["Terpene"] = (0.2, 0.75, 0.05, 2.0)
    bgc_class_weight["PKS-NRP_Hybrids"] = (0.0, 0.78, 0.22, 1.0)
    bgc_class_weight["Others"] = (0.01, 0.97, 0.02, 4.0)
    
    #define which classes will be analyzed (if in the options_classify mode)
    valid_classes = set()
    for key in bgc_class_weight:
        valid_classes.add(key.lower())
    user_banned_classes = set([a.strip().lower() for a in options.banned_classes])
    valid_classes = valid_classes - user_banned_classes
    
    # finally, define weights for mix
    bgc_class_weight["mix"] = (0.2, 0.75, 0.05, 2.0) # default when not separating in classes
    BGC_classes = defaultdict(list)
    # mix class will always be the last element of the tuple
    bgcClassNames = tuple(sorted(list(bgc_class_weight)) + ["mix"])
    assert bgcClassNames[-1] == 'mix'


    # genbankDict: {cluster_name:[genbank_path_to_1st_instance,[sample_1,sample_2,...]]}
    bgc_info = {} # Stores, per BGC: predicted type, gbk Description, number of records, width of longest record, GenBank's accession, Biosynthetic Genes' ids
    genbankDict = {}
    
    # Exclude single string
    include_gbk_str = options.include_gbk_str
    if len(include_gbk_str) == 1 and include_gbk_str[0] == "*":
        print(" Including all files")
    elif len(include_gbk_str) == 1 and include_gbk_str[0] == "":
        sys.exit(" Stop: no strings specified for '--include_gbk_str'")
    else:
        print(" Including files with one or more of the following strings in their filename: '{}'".format("', '".join(include_gbk_str)))
    
    exclude_gbk_str = options.exclude_gbk_str
    if exclude_gbk_str != []:
        print(" Skipping files with one or more of the following strings in their filename: '{}'".format("', '".join(exclude_gbk_str)))
    
    # Read included MIBiG
    # Change this for every officially curated MIBiG bundle
    # (file, final folder, number of bgcs)
    mibig_set = set()
    if use_relevant_mibig:
        if options.mibig21:
            mibig_zipfile_numbgcs = ("MIBiG_2.1_final.zip", "MIBiG_2.1_final", 1923)
        elif options.mibig14:
            mibig_zipfile_numbgcs = ("MIBiG_1.4_final.zip", "MIBiG_1.4_final", 1808)
        else:
            mibig_zipfile_numbgcs = ("MIBiG_1.3_final.zip", "MIBiG_1.3_final", 1393)
        
        print("\n Trying to read bundled MIBiG BGCs as reference")
        mibig_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),"Annotated_MIBiG_reference")
        bgcs_path = os.path.join(mibig_path,mibig_zipfile_numbgcs[1])
        
        # try to see if the zip file has already been decompressed
        numbgcs = len(glob(os.path.join(bgcs_path,"*.gbk")))
        if numbgcs == 0:
            if not zipfile.is_zipfile(os.path.join(mibig_path,mibig_zipfile_numbgcs[0])):
                sys.exit("Did not find file {}. Please re-download it from the official repository".format(mibig_zipfile_numbgcs[0]))
                
            with zipfile.ZipFile(os.path.join(mibig_path,mibig_zipfile_numbgcs[0]), 'r') as mibig_zip:
                for fname in mibig_zip.namelist():
                    if fname[-3:] != "gbk":
                        continue
                
                    extractedbgc = mibig_zip.extract(fname,path=mibig_path)
                    if verbose:
                        print("  Extracted {}".format(extractedbgc))
        
        elif mibig_zipfile_numbgcs[2] == numbgcs:
            print("  MIBiG BGCs seem to have been extracted already")
        else:
            sys.exit("Did not find the correct number of MIBiG BGCs ({}). Please clean the 'Annotated MIBiG reference' folder from any .gbk files first".format(mibig_zipfile_numbgcs[2]))
        
        print("\nImporting MIBiG files")
        get_gbk_files(bgcs_path, output_folder, bgc_fasta_folder, int(options.min_bgc_size),
                      ['*'], exclude_gbk_str, bgc_info, mode, verbose, force_hmmscan, valid_classes, bgc_data, genbankDict)
        
        for i in genbankDict.keys():
            mibig_set.add(i)
            
    
    print("\nImporting GenBank files")
    get_gbk_files(options.inputdir, output_folder, bgc_fasta_folder, int(options.min_bgc_size),
                  include_gbk_str, exclude_gbk_str, bgc_info, mode, verbose, force_hmmscan, valid_classes, bgc_data, genbankDict)
    
    if has_query_bgc:
        query_bgc = ".".join(options.query_bgc.split(os.sep)[-1].split(".")[:-1])
        if query_bgc in genbankDict:
            print("\nQuery BGC already added")
            pass
        else:
            print("\nImporting query BGC file")
            get_gbk_files(options.query_bgc, output_folder, bgc_fasta_folder, 
                          int(options.min_bgc_size), ['*'], exclude_gbk_str, bgc_info)
            
        if query_bgc not in genbankDict:
            sys.exit("Error: not able to include Query BGC (check valid classes, BGC size, etc. Run again with --verbose)")
    # clusters and sampleDict contain the necessary structure for all-vs-all and sample analysis
    clusters = list(genbankDict.keys())
    
    sampleDict = {} # {sampleName:set(bgc1,bgc2,...)}
    gbk_files = [] # raw list of gbk file locations
    for (cluster, (path, clusterSample)) in genbankDict.items():
        gbk_files.append(path)
        for sample in clusterSample:
            clustersInSample = sampleDict.get(sample, set())
            clustersInSample.add(cluster)
            sampleDict[sample] = clustersInSample
    
    print("\nCreating output directories")
    svg_folder = os.path.join(output_folder, "SVG")
    create_directory(svg_folder, "SVG", False)
    network_folder = os.path.join(output_folder, "network_files")
    create_directory(network_folder, "Networks", False)

    print("\nTrying threading on {} cores".format(str(cores)))
    
    """BGCs -- 
    dictionary of this structure:
    BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
     'specific_domain_name_2'] } }
    - cluster_name_x: cluster name (can be anything)
    - general_domain_name_x: PFAM ID, for example 'PF00550'
    - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names in DMS unequivocally
     (for example, 'PF00550_start_end', where start and end are genomic positions)."""
    BGCs = {} #will contain the BGCs

    bgcClassName2idx = dict(zip(bgcClassNames,range(len(bgcClassNames))))

    AlignedDomainSequences = {} # Key: specific domain sequence label. Item: aligned sequence
    DomainList = {} # Key: BGC. Item: ordered list of domains
    
    # Key: BGC. Item: ordered list of simple integers with the number of domains
    # in each gene
    # Instead of `DomainCountGene = defaultdict(list)`, let's try arrays of 
    # unsigned ints
    DomainCountGene = {}
    # list of gene-numbers that have a hit in the anchor domain list. Zero based
    corebiosynthetic_position = {}
    # list of +/- orientation 
    BGCGeneOrientation = {}
    
    # to avoid multiple alignment if there's only 1 seq. representing a particular domain
    sequences_per_domain = {}
    
    
    ### Step 2: Run hmmscan
    print("\nPredicting domains using hmmscan")
    
    baseNames = set(clusters)
    
    # All available fasta files (could be more than it should if reusing output folder)
    allFastaFiles = set(glob(os.path.join(bgc_fasta_folder,"*.fasta")))
    
    # fastaFiles: all the fasta files that should be there 
    # (i.e. correspond to the input files)
    fastaFiles = set()
    for name in baseNames:
        fastaFiles.add(os.path.join(bgc_fasta_folder, name+".fasta"))

    # fastaBases: the actual fasta files we have that correspond to the input
    fastaBases = allFastaFiles.intersection(fastaFiles)
    
    # Verify that all input files had their fasta sequences extracted
    if len(fastaFiles - fastaBases) > 0:
        sys.exit("Error! The following files did NOT have their fasta sequences extracted: " + ", ".join(fastaFiles - fastaBases))
    
    # Make a list of all fasta files that need to be processed
    # (i.e., they don't yet have a corresponding .domtable)
    if force_hmmscan:
        # process all files, regardless of whether they already existed
        task_set = fastaFiles
        print(" Forcing domain prediction on ALL fasta files (--force_hmmscan)")
    else:
        # find already processed files
        alreadyDone = set()
        for fasta in fastaFiles:
            outputbase  = ".".join(fasta.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(domtable_folder,outputbase + '.domtable')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                # verify domtable content
                with open(outputfile, "r") as domtablefile:
                    for line in domtablefile.readlines():
                        if line.startswith("# Option settings:"):
                            linecols = line.split()
                            if "hmmscan" in linecols and "--domtblout" in linecols:
                                alreadyDone.add(fasta)
                                break
                
        task_set = fastaFiles - alreadyDone
        if len(task_set) == 0:
            print(" All fasta files had already been processed")
        elif len(alreadyDone) > 0:
            if len(task_set) < 20:
                print(" Warning! The following NEW fasta file(s) will be processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in task_set)))
            else:
                print(" Warning: {} NEW fasta files will be processed".format(str(len(task_set))))
        else:
            print(" Predicting domains for {} fasta files".format(str(len(fastaFiles))))
        
    pool = Pool(cores,maxtasksperchild=1)
    for fastaFile in task_set:
        pool.apply_async(runHmmScan,args=(fastaFile, pfam_dir, domtable_folder, verbose))
    pool.close()
    pool.join()
    print(" Finished generating domtable files.")

    ### Step 3: Parse hmmscan domtable results and generate pfs and pfd files
    print("\nParsing hmmscan domtable files")
    
    # All available domtable files
    allDomtableFiles = set(glob(os.path.join(domtable_folder,"*.domtable")))
    
    # domtableFiles: all domtable files corresponding to the input files
    domtableFiles = set()
    for name in baseNames:
        domtableFiles.add(os.path.join(domtable_folder, name+".domtable"))
    
    # domtableBases: the actual set of input files with coresponding domtable files
    domtableBases = allDomtableFiles.intersection(domtableFiles)
    
    # Verify that all input files have a corresponding domtable file
    if len(domtableFiles - domtableBases) > 0:
        sys.exit("Error! The following files did NOT have their domains predicted: " + ", ".join(domtableFiles - domtableBases))
    
    # find already processed files (assuming that if the pfd file exists, the pfs should too)
    alreadyDone = set()
    if not force_hmmscan:
        for domtable in domtableFiles:
            outputbase = ".".join(domtable.split(os.sep)[-1].split(".")[:-1])
            outputfile = os.path.join(pfd_folder, outputbase + '.pfd')
            if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0:
                alreadyDone.add(domtable)
    domtableFilesUnprocessed = domtableFiles - alreadyDone
    if len(domtableFilesUnprocessed) == 0: # Re-run
        print(" All domtable files had already been processed")
    elif len(alreadyDone) > 0: # Incomplete run
        if len(domtableFilesUnprocessed) < 20:
            print(" Warning! The following domtable files had not been processed: {}".format(", ".join(".".join(x.split(os.sep)[-1].split(".")[:-1]) for x in domtableFilesUnprocessed)))
        else:
            print(" Warning: {} domtable files will be processed".format(str(len(domtableFilesUnprocessed))))
    else: # First run
        print(" Processing {} domtable files".format(str(len(domtableFiles))))

    # If using the multiprocessing version and outputbase doesn't have any
    #  predicted domains, it's not as easy to remove if from the analysis
    #  (probably because parseHmmScan only has a copy of clusters et al?)
    # Using serialized version for now. Probably doesn't have too bad an impact
    #pool = Pool(cores,maxtasksperchild=32)
    for domtableFile in domtableFiles - alreadyDone:
        parseHmmScan(domtableFile, pfd_folder, pfs_folder, options.domain_overlap_cutoff, verbose, genbankDict, clusters, baseNames, gbk_files, sampleDict, mibig_set)
        #pool.apply_async(parseHmmScan, args=(domtableFile,output_folder,options.domain_overlap_cutoff))
    #pool.close()
    #pool.join()
    
    # If number of pfd files did not change, no new sequences were added to the 
    #  domain fastas and we could try to resume the multiple alignment phase
    # baseNames have been pruned of BGCs with no domains that might've been added temporarily
    try_MA_resume = False
    if len(baseNames - set(pfd.split(os.sep)[-1][:-9] for pfd in alreadyDone)) == 0:
        try_MA_resume = True
    else:
        # new sequences will be added to the domain fasta files. Clean domains folder
        # We could try to make it so it's not necessary to re-calculate every alignment,
        #  either by expanding previous alignment files or at the very least, 
        #  re-aligning only the domain files of the newly added BGCs
        print(" New domain sequences to be added; cleaning domains folder")
        for thing in os.listdir(domains_folder):
            os.remove(os.path.join(domains_folder,thing))

    print(" Finished generating pfs and pfd files.")
    

    ### Step 4: Parse the pfs, pfd files to generate BGC dictionary, clusters, and clusters per sample objects
    print("\nProcessing domains sequence files")
    
    # All available pfd files
    allPfdFiles = set(glob(os.path.join(pfd_folder,"*.pfd")))
    
    # pfdFiles: all pfd files corresponding to the input files
    # (some input files could've been removed due to not having predicted domains)
    pfdFiles = set()
    for name in baseNames:
        pfdFiles.add(os.path.join(pfd_folder, name+".pfd"))
    
    # pfdBases: the actual set of input files that have pfd files
    pfdBases = allPfdFiles.intersection(pfdFiles)
    
    # verify previous step. 
    # All BGCs without predicted domains should no longer be in baseNames
    if len(pfdFiles - pfdBases) > 0:
        sys.exit("Error! The following files did NOT have their domtable files processed: " + ", ".join(pfdFiles - pfdBases))

    filtered_matrix = []
    if options.skip_ma:
        print(" Running with skip_ma parameter: Assuming that the domains folder has all the fasta files")
        try:
            with open(os.path.join(cache_folder, "BGCs.dict"), "r") as BGC_file:
                BGCs = pickle.load(BGC_file)
                BGC_file.close()
        except IOError:
            sys.exit("BGCs file not found...")
    else:
        print(" Adding sequences to corresponding domains file")
            
        for outputbase in baseNames:
            if verbose:
                print("   Processing: " + outputbase)

            pfdFile = os.path.join(pfd_folder, outputbase + ".pfd")
            filtered_matrix = [[part.strip() for part in line.split('\t')] for line in open(pfdFile)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(bgc_fasta_folder, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and 
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not try_MA_resume:
                with open(fasta_file, "r") as fasta_file_handle:
                    fasta_dict = fasta_parser(fasta_file_handle) # all fasta info from a BGC
                save_domain_seqs(filtered_matrix, fasta_dict, domains_folder, outputbase)

            BGCs[outputbase] = BGC_dic_gen(filtered_matrix)
            
            del filtered_matrix[:]
            
        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(cache_folder, "BGCs.dict"), "wb") as BGC_file:
            pickle.dump(BGCs, BGC_file)
            BGC_file.close()
            
    # if it's a re-run, the pfd/pfs files were not changed, so the skip_ma flag
    # is activated. We have to open the pfd files to get the gene labels for
    # each domain
    # We now always have to have this data so the alignments are produced
    pfd_dict_domains = defaultdict(int)
    orf_keys = {}
    for outputbase in baseNames:
        DomainCountGene[outputbase] = array('B')
        corebiosynthetic_position[outputbase] = array('H')
        BGCGeneOrientation[outputbase] = array('b')
        pfdFile = os.path.join(pfd_folder, outputbase + ".pfd")
        
        #pfd_dict_domains contains the number of domains annotated in the
        # pfd file for each orf tag
        with open(pfdFile,"r") as pfdf:
            for line in pfdf:
                pfd_dict_domains[line.strip().split("\t")[-1]] += 1
        
        # extract the orf number from the tag and use it to traverse the BGC
        for orf in pfd_dict_domains.keys():
            orf_num = int(orf.split(":")[0].split("_ORF")[1])
            orf_keys[orf_num] = orf
        
        orf_num = 0
        for orf_key in sorted(orf_keys.keys()):
            orf = orf_keys[orf_key]
            if orf[-1] == "+":
                BGCGeneOrientation[outputbase].append(1)
            else:
                BGCGeneOrientation[outputbase].append(-1)
                
            DomainCountGene[outputbase].append(pfd_dict_domains[orf])
            
            if orf in bgc_info[outputbase].biosynthetic_genes:
                corebiosynthetic_position[outputbase].append(orf_num)
            orf_num += 1
        
        pfd_dict_domains.clear()
        orf_keys.clear()

        ## TODO: if len(corebiosynthetic_position[outputbase]) == 0
        ## do something with the list of pfam ids. Specifically, mark
        ## (in this case TODO or always?) as biosynthetic genes, the ones that contain
        ## domains from a special list. This list of special domains
        ## comes from predicted domains within the CDSs marked as 'sec_met'
        ## by antismash
            
            
    # Get the ordered list of domains
    print(" Reading the ordered list of domains from the pfs files")
    for outputbase in baseNames:
        pfsfile = os.path.join(pfs_folder, outputbase + ".pfs")
        if os.path.isfile(pfsfile):
            DomainList[outputbase] = get_domain_list(pfsfile)
        else:
            sys.exit(" Error: could not open " + outputbase + ".pfs")
                
                
    ### Step 5: Create SVG figures
    print(" Creating arrower-like figures for each BGC")
    
    # read hmm file. We'll need that info anyway for final visualization
    print("  Parsing hmm file for domain information")
    pfam_info = {}
    with open(os.path.join(pfam_dir, "Pfam-A.hmm"), "r") as pfam:
        putindict = False
        # assuming that the order of the information never changes
        for line in pfam:
            if line[:4] == "NAME":
                name = line.strip()[6:]
            if line[:3] == "ACC":
                acc = line.strip()[6:].split(".")[0]
            if line[:4] == "DESC":
                desc = line.strip()[6:]
                putindict = True
                
            if putindict:
                putindict = False
                pfam_info[acc] = (name, desc)
    print("    Done")
    
    # verify if there are figures already generated
    
    # All available SVG files
    availableSVGs = set()
    for svg in glob(os.path.join(svg_folder,"*.svg")):
        (root, ext) = os.path.splitext(svg)
        availableSVGs.add(root.split(os.sep)[-1])
        
    # Which files actually need to be generated
    working_set = baseNames - availableSVGs
    
    if len(working_set) > 0:
        color_genes = {}
        color_domains = read_color_domains_file()
        pfam_domain_categories = {}
        
        #This must be done serially, because if a color for a gene/domain
        # is not found, the text files with colors need to be updated
        print("  Reading BGC information and writing SVG")
        for bgc in working_set:
            with open(genbankDict[bgc][0],"r") as handle:
                SVG(False, os.path.join(svg_folder,bgc+".svg"), handle, bgc, os.path.join(pfd_folder,bgc+".pfd"), True, color_genes, color_domains, pfam_domain_categories, pfam_info, bgc_info[bgc].records, bgc_info[bgc].max_width)
        
        color_genes.clear()
        color_domains.clear()
        pfam_domain_categories.clear()
    elif len(working_set) == 0:
        print("  All SVG from the input files seem to be in the SVG folder")
    
    availableSVGs.clear()
    print(" Finished creating figures")
    
    
    print("\n\n   - - Calculating distance matrix - -")
   
    # Do multiple alignments if needed
    if not options.skip_ma:
        print("Performing multiple alignment of domain sequences")
        
        # obtain all fasta files with domain sequences
        domain_sequence_list = set(glob(os.path.join(domains_folder,"*.fasta")))
        
        # compare with .algn set of files. Maybe resuming is possible if
        # no new sequences were added
        if try_MA_resume:
            temp_aligned = set(glob(os.path.join(domains_folder, "*.algn")))
            
            if len(temp_aligned) > 0:
                print(" Found domain fasta files without corresponding alignments")
                
                for a in temp_aligned:
                    if os.path.getsize(a) > 0:
                        domain_sequence_list.remove(a[:-5]+".fasta")
            
            temp_aligned.clear()
        
        # Try to further reduce the set of domain fastas that need alignment
        sequence_tag_list = set()
        header_list = []
        domain_sequence_list_temp = domain_sequence_list.copy()
        for domain_file in domain_sequence_list_temp:
            domain_name = ".".join(domain_file.split(os.sep)[-1].split(".")[:-1])
            
            # fill fasta_dict...
            with open(domain_file, "r") as fasta_handle:
                header_list = get_fasta_keys(fasta_handle)
                
            # Get the BGC name from the sequence tag. The form of the tag is:
            # >BGCXXXXXXX_BGCXXXXXXX_ORF25:gid...
            sequence_tag_list = set(s.split("_ORF")[0] for s in header_list)

            # ...to find out how many sequences do we actually have
            if len(sequence_tag_list) == 1:
                # avoid multiple alignment if the domains all belong to the same BGC
                domain_sequence_list.remove(domain_file)
                if verbose:
                    print(" Skipping Multiple Alignment for {} (appears only in one BGC)".format(domain_name))
        
        sequence_tag_list.clear()
        del header_list[:]
        
        domain_sequence_list_temp.clear()
            
        # Do the multiple alignment
        stop_flag = False
        if len(domain_sequence_list) > 0:
            print("\n Using hmmalign")
            launch_hmmalign(cores, domain_sequence_list, pfam_dir, verbose)
                
            # verify all tasks were completed by checking existance of alignment files
            for domain_file in domain_sequence_list:
                if not os.path.isfile(domain_file[:-6]+".algn"):
                    print("   ERROR, {}.algn could not be found (possible issue with aligner).".format(domain_file[:-6]))
                    stop_flag = True
            if stop_flag:
                sys.exit()
                       
        else:
            print(" No domain fasta files found to align")
    
    
    # If there's something to analyze, load the aligned sequences
    print(" Trying to read domain alignments (*.algn files)")
    aligned_files_list = glob(os.path.join(domains_folder, "*.algn"))
    if len(aligned_files_list) == 0:
        sys.exit("No aligned sequences found in the domain folder (run without the --skip_ma parameter or point to the correct output folder)")
    for aligned_file in aligned_files_list:
        with open(aligned_file, "r") as aligned_file_handle:
            fasta_dict = fasta_parser(aligned_file_handle)
            for header in fasta_dict:
                AlignedDomainSequences[header] = fasta_dict[header]

    clusterNames = tuple(sorted(clusters))
    
    # we have to find the idx of query_bgc
    if has_query_bgc:
        try:
            query_bgc_idx = clusterNames.index(query_bgc)
        except ValueError:
            sys.exit("Error finding the index of Query BGC")

    # create output directory for network files
    network_files_folder = os.path.join(network_folder, run_name)
    create_directory(network_files_folder, "Network Files", False)

    # copy html templates
    dir_util.copy_tree(os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "output"), output_folder)

    # make a new run folder in the html output & copy the overview_html
    network_html_folder = os.path.join(output_folder, "html_content", "networks", run_name)
    rundata_networks_per_run = {}
    html_subs_per_run = {}
    for cutoff in cutoff_list:
        network_html_folder_cutoff = "{}_c{:.2f}".format(network_html_folder, cutoff)
        create_directory(network_html_folder_cutoff, "Network HTML Files", False)
        shutil.copy(os.path.join(os.path.dirname(os.path.realpath(__file__)), "html_template", "overview_html"), os.path.join(network_html_folder_cutoff, "overview.html"))
        rundata_networks_per_run[network_html_folder_cutoff] = []
        html_subs_per_run[network_html_folder_cutoff] = []
        
    # create pfams.js
    pfams_js_file = os.path.join(output_folder, "html_content", "js", "pfams.js")
    if not os.path.isfile(pfams_js_file):
        with open(pfams_js_file, "w") as pfams_js:
            pfam_json = {}
            pfam_colors = generatePfamColorsMatrix(os.path.join(os.path.dirname(os.path.realpath(__file__)), "domains_color_file.tsv"))
            for pfam_code in pfam_info:
                pfam_obj = {}
                if pfam_code in pfam_colors:
                    pfam_obj["col"] = pfam_colors[pfam_code]
                else:
                    pfam_obj["col"] = "255,255,255"
                pfam_obj["desc"] = pfam_info[pfam_code][1]
                pfam_json[pfam_code] = pfam_obj
            pfams_js.write("var pfams={};\n".format(json.dumps(pfam_json, indent=4, separators=(',', ':'), sort_keys=True)))

    # Try to make default analysis using all files found inside the input folder
    print("\nGenerating distance network files with ALL available input files")

    # This version contains info on all bgcs with valid classes
    print("   Writing the complete Annotations file for the complete set")
    network_annotation_path = os.path.join(network_files_folder, "Network_Annotations_Full.tsv")
    with open(network_annotation_path, "w") as network_annotation_file:
        network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
        for bgc in clusterNames:
            product = bgc_info[bgc].product
            network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
    
    
    # Find index of all MIBiG BGCs if necessary
    if use_relevant_mibig:
        name_to_idx = {}
        for clusterIdx,clusterName in enumerate(clusterNames):
            name_to_idx[clusterName] = clusterIdx
            
        mibig_set_indices = set()
        for bgc in mibig_set:
            mibig_set_indices.add(name_to_idx[bgc])

    # Making network files mixing all classes
    if options_mix:
        print("\n Mixing all BGC classes")
        
        # only choose from valid classes
        mix_set = []
        
        # create working set with indices of valid clusters
        for clusterIdx,clusterName in enumerate(clusterNames):
            if has_includelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DomainList[clusterName]})
                    
                if len(domain_includelist & bgc_domain_set) == 0:
                    continue
            
            product = bgc_info[clusterName].product
            predicted_class = sort_bgc(product)
            
            if predicted_class.lower() in valid_classes:
                mix_set.append(clusterIdx)
        
        print("\n  {} ({} BGCs)".format("Mix", str(len(mix_set))))

        # create output directory
        create_directory(os.path.join(network_files_folder, "mix"), "  Mix", False)
        
        print("  Calculating all pairwise distances")
        if has_query_bgc:
            pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx], mix_set)])
        else:
            # convert into a set of ordered tuples
            pairs = set([tuple(sorted(combo)) for combo in combinations(mix_set, 2)])
        
        cluster_pairs = [(x, y, -1) for (x, y) in pairs]
        pairs.clear()
        network_matrix_mix = generate_network(cluster_pairs, cores, clusterNames, bgcClassNames, DomainList, output_folder, DomainCountGene, 
        corebiosynthetic_position, BGCGeneOrientation, bgc_class_weight, anchor_domains, BGCs, mode, bgc_info,
        AlignedDomainSequences, verbose, domains_folder)
        
        del cluster_pairs[:]

        # add the rest of the edges in the "Query network"
        if has_query_bgc:
            new_set = []
            
            # rows from the distance matrix that will be pruned 
            del_list = [] 
        
            for idx, row in enumerate(network_matrix_mix):
                a, b, distance = int(row[0]), int(row[1]), row[2]
                
                if a == b:
                    continue
                
                if distance <= max_cutoff:
                    if a == query_bgc_idx:
                        new_set.append(b)
                    else:
                        new_set.append(a)
                else:
                    del_list.append(idx)
            
            for idx in sorted(del_list, reverse=True):
                del network_matrix_mix[idx]
            del del_list[:]
            
            pairs = set([tuple(sorted(combo)) for combo in combinations(new_set, 2)])
            cluster_pairs = [(x, y, -1) for (x, y) in pairs]
            pairs.clear()
            network_matrix_new_set = generate_network(cluster_pairs, cores, clusterNames, bgcClassNames, DomainList, output_folder, DomainCountGene,
            corebiosynthetic_position, BGCGeneOrientation, bgc_class_weight, anchor_domains, BGCs, mode, bgc_info,
            AlignedDomainSequences, verbose, domains_folder)
            del cluster_pairs[:]
            
            # Update the network matrix (QBGC-vs-all) with the distances of
            # QBGC's GCF
            network_matrix_mix.extend(network_matrix_new_set)
            
            # Update actual list of BGCs that we'll use
            mix_set = new_set
            mix_set.extend([query_bgc_idx])
            mix_set.sort() # clusterJsonBatch expects ordered indices
        
            # Create an additional file with the list of all clusters in the class + other info
            # This version of the file only has information on the BGCs connected to Query BGC
            print("   Writing annotation file")
            network_annotation_path = os.path.join(network_files_folder, "mix", "Network_Annotations_mix_QueryBGC.tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in mix_set:
                    bgc = clusterNames[idx]
                    product = bgc_info[bgc].product
                    network_annotation_file.write("\t".join([bgc, 
                        bgc_info[bgc].accession_id, bgc_info[bgc].description, 
                        product, sort_bgc(product), bgc_info[bgc].organism, 
                        bgc_info[bgc].taxonomy]) + "\n")
        elif use_relevant_mibig:
            n = nx.Graph()
            n.add_nodes_from(mix_set)
            mibig_set_del = []
            network_matrix_set_del = []
            
            for idx, row in enumerate(network_matrix_mix):
                a, b, distance = int(row[0]), int(row[1]), row[2]
                if distance <= max_cutoff:
                    n.add_edge(a, b, index=idx)
                    
            for component in nx.connected_components(n): # note: 'component' is a set
                numBGCs_subgraph = len(component)
                
                # catch if the subnetwork is comprised only of MIBiG BGCs
                if len(component & mibig_set_indices) == numBGCs_subgraph:
                    for bgc in component:
                        mibig_set_del.append(bgc)
                    
            # Get all edges between bgcs marked for deletion
            for (a, b, idx) in n.subgraph(mibig_set_del).edges.data('index'):
                network_matrix_set_del.append(idx)
                
            # delete all edges between marked bgcs
            for row_idx in sorted(network_matrix_set_del, reverse=True):
                del network_matrix_mix[row_idx]
            del network_matrix_set_del[:]
            
            print("   Removing {} non-relevant MIBiG BGCs".format(len(mibig_set_del)))
            mix_set_idx = 0
            bgc_to_mix_set_idx = {}
            for idx, bgc in enumerate(mix_set):
                bgc_to_mix_set_idx[bgc] = idx
            
            for bgc_idx in sorted(mibig_set_del, reverse=True):
                del mix_set[bgc_to_mix_set_idx[bgc_idx]]
            del mibig_set_del[:]
            

        print("  Writing output files")
        pathBase = os.path.join(network_files_folder, "mix")
        filenames = []
        for cutoff in cutoff_list:
            filenames.append(os.path.join(pathBase, "mix_c{:.2f}.network".format(cutoff)))
        cutoffs_and_filenames = list(zip(cutoff_list, filenames))
        del filenames[:]
        write_network_matrix(network_matrix_mix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)
            
        print("  Calling Gene Cluster Families")
        reduced_network = []
        pos_alignments = {}
        for row in network_matrix_mix:
            reduced_network.append([int(row[0]), int(row[1]), row[2]])
            reverse = False
            if row[-1] == 1.0:
                reverse = True
            pa = pos_alignments.setdefault(int(row[0]),{})
            # lcsStartA, lcsStartB, seedLength, reverse={True,False}
            pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
        del network_matrix_mix[:]
        family_data = clusterJsonBatch(mix_set, pathBase, "mix", reduced_network, pos_alignments,
            clusterNames, bgc_info, mibig_set, pfd_folder, bgc_fasta_folder,
            DomainList, BGCs, AlignedDomainSequences, DomainCountGene, BGCGeneOrientation,
            cutoffs=cutoff_list, clusterClans=options.clans,
            clanCutoff=options.clan_cutoff, htmlFolder=network_html_folder)
        for network_html_folder_cutoff in family_data:
            rundata_networks_per_run[network_html_folder_cutoff].append(family_data[network_html_folder_cutoff])
            html_subs_per_run[network_html_folder_cutoff].append({ "name" : "mix", "css" : "Others", "label" : "Mixed"})
        del mix_set[:]
        del reduced_network[:]
        
        
    # Making network files separating by BGC class
    if options_classify:
        print("\n Working for each BGC class")
        
        # reinitialize BGC_classes to make sure the bgc lists are empty
        BGC_classes = defaultdict(list)
    
        # Preparing gene cluster classes
        print("  Sorting the input BGCs\n")
        
        # create and sort working set for each class
        for clusterIdx,clusterName in enumerate(clusterNames):
            if has_includelist:
                # extra processing because pfs info includes model version
                bgc_domain_set = set({x.split(".")[0] for x in DomainList[clusterName]})
                    
                if len(domain_includelist & bgc_domain_set) == 0:
                    continue
            
            product = bgc_info[clusterName].product
            predicted_class = sort_bgc(product)
            
            if predicted_class.lower() in valid_classes:
                BGC_classes[predicted_class].append(clusterIdx)
            
            # possibly add hybrids to 'pure' classes
            if options.hybrids:
                if predicted_class == "PKS-NRP_Hybrids":
                    if "nrps" in valid_classes:
                        BGC_classes["NRPS"].append(clusterIdx)
                    if "t1pks" in product and "pksi" in valid_classes:
                        BGC_classes["PKSI"].append(clusterIdx)
                    if "t1pks" not in product and "pksother" in valid_classes:
                        BGC_classes["PKSother"].append(clusterIdx)
                
                if predicted_class == "Others" and "." in product:
                    subclasses = set()
                    for subproduct in product.split("."):
                        subclass = sort_bgc(subproduct)
                        if subclass.lower() in valid_classes:
                            subclasses.add(subclass)
                            
                    # Prevent mixed BGCs with sub-Others annotations to get
                    # added twice (e.g. indole-cf_fatty_acid has already gone
                    # to Others at this point)
                    if "Others" in subclasses:
                        subclasses.remove("Others")
                        
                        
                    for subclass in subclasses:
                        BGC_classes[subclass].append(clusterIdx)
                    subclasses.clear()

        # only make folders for the BGC_classes that are found
        for bgc_class in BGC_classes:
            if has_query_bgc:
                # not interested in this class if our Query BGC is not here...
                if query_bgc_idx not in BGC_classes[bgc_class]:
                    continue
            
            print("\n  {} ({} BGCs)".format(bgc_class, str(len(BGC_classes[bgc_class]))))
            if use_relevant_mibig:
                if len(set(BGC_classes[bgc_class]) & mibig_set_indices) == len(BGC_classes[bgc_class]):
                    print(" - All clusters in this class are MIBiG clusters -")
                    print("  If you'd like to analyze MIBiG clusters, turn off the --mibig option")
                    print("  and point --inputdir to the Annotated_MIBiG_reference folder")
                    continue
            
            # create output directory
            create_directory(os.path.join(network_files_folder, bgc_class), "  All - " + bgc_class, False)
            
            # Create an additional file with the final list of all clusters in the class
            print("   Writing annotation files")
            network_annotation_path = os.path.join(network_files_folder, bgc_class, "Network_Annotations_" + bgc_class + ".tsv")
            with open(network_annotation_path, "w") as network_annotation_file:
                network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                for idx in BGC_classes[bgc_class]:
                    bgc = clusterNames[idx]
                    product = bgc_info[bgc].product
                    network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
            
            print("   Calculating all pairwise distances")
            if has_query_bgc:
                pairs = set([tuple(sorted(combo)) for combo in combinations_product([query_bgc_idx],BGC_classes[bgc_class])])
            else:
                pairs = set([tuple(sorted(combo)) for combo in combinations(BGC_classes[bgc_class], 2)])
                
            cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
            pairs.clear()
            network_matrix = generate_network(cluster_pairs, cores, clusterNames, bgcClassNames, DomainList, output_folder, DomainCountGene,
            corebiosynthetic_position, BGCGeneOrientation, bgc_class_weight, anchor_domains, BGCs, mode, bgc_info,
            AlignedDomainSequences, verbose, domains_folder)
            #pickle.dump(network_matrix,open("others.ntwrk",'wb'))
            del cluster_pairs[:]
            #network_matrix = pickle.load(open("others.ntwrk", "rb"))
                
            # add the rest of the edges in the "Query network"
            if has_query_bgc:
                new_set = []
                
                # rows from the distance matrix that will be pruned 
                del_list = []
                
                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]
                    
                    # avoid QBGC-QBGC
                    if a == b:
                        continue
                    
                    if distance <= max_cutoff:
                        if a == query_bgc_idx:
                            new_set.append(b)
                        else:
                            new_set.append(a)
                    else:
                        del_list.append(idx)
                
                for idx in sorted(del_list, reverse=True):
                    del network_matrix[idx]
                del del_list[:]
                
                pairs = set([tuple(sorted(combo)) for combo in combinations(new_set, 2)])
                cluster_pairs = [(x, y, bgcClassName2idx[bgc_class]) for (x, y) in pairs]
                pairs.clear()
                network_matrix_new_set = generate_network(cluster_pairs, cores, clusterNames, bgcClassNames, DomainList, output_folder, DomainCountGene,
                corebiosynthetic_position, BGCGeneOrientation, bgc_class_weight, anchor_domains, BGCs, mode, bgc_info,
                AlignedDomainSequences, verbose, domains_folder)
                del cluster_pairs[:]
                                    
                # Update the network matrix (QBGC-vs-all) with the distances of
                # QBGC's GCF
                network_matrix.extend(network_matrix_new_set)
                
                # Update actual list of BGCs that we'll use
                BGC_classes[bgc_class] = new_set
                BGC_classes[bgc_class].extend([query_bgc_idx])
                BGC_classes[bgc_class].sort()
                
                # Create an additional file with the list of all clusters in the class + other info
                # This version of the file only has information on the BGCs connected to Query BGC
                print("   Writing annotation file (Query BGC)")
                network_annotation_path = os.path.join(network_files_folder, bgc_class, "Network_Annotations_" + bgc_class + "_QueryBGC.tsv")
                with open(network_annotation_path, "w") as network_annotation_file:
                    network_annotation_file.write("BGC\tAccession ID\tDescription\tProduct Prediction\tBiG-SCAPE class\tOrganism\tTaxonomy\n")
                    for idx in BGC_classes[bgc_class]:
                        bgc = clusterNames[idx]
                        product = bgc_info[bgc].product
                        network_annotation_file.write("\t".join([bgc, bgc_info[bgc].accession_id, bgc_info[bgc].description, product, sort_bgc(product), bgc_info[bgc].organism, bgc_info[bgc].taxonomy]) + "\n")
            elif use_relevant_mibig:
                n = nx.Graph()
                n.add_nodes_from(BGC_classes[bgc_class])
                mibig_set_del = []
                network_matrix_set_del = []
                
                for idx, row in enumerate(network_matrix):
                    a, b, distance = int(row[0]), int(row[1]), row[2]
                    if distance <= max_cutoff:
                        n.add_edge(a, b, index=idx)
                        
                for component in nx.connected_components(n): # note: 'component' is a set
                    numBGCs_subgraph = len(component)
                    
                    # catch if the subnetwork is comprised only of MIBiG BGCs
                    if len(component & mibig_set_indices) == numBGCs_subgraph:
                        for bgc in component:
                            mibig_set_del.append(bgc)
                
                # Get all edges between bgcs marked for deletion
                for (a, b, idx) in n.subgraph(mibig_set_del).edges.data('index'):
                    network_matrix_set_del.append(idx)
                
                # delete all edges between marked bgcs
                for row_idx in sorted(network_matrix_set_del, reverse=True):
                    del network_matrix[row_idx]
                del network_matrix_set_del[:]
                            
                print("   Removing {} non-relevant MIBiG BGCs".format(len(mibig_set_del)))
                bgc_to_class_idx = {}
                for idx, bgc in enumerate(BGC_classes[bgc_class]):
                    bgc_to_class_idx[bgc] = idx
                for bgc_idx in sorted(mibig_set_del, reverse=True):
                    del BGC_classes[bgc_class][bgc_to_class_idx[bgc_idx]]
                del mibig_set_del[:]
                    
                
                
            if len(BGC_classes[bgc_class]) < 2:
                continue
                
            print("   Writing output files")
            pathBase = os.path.join(network_files_folder, bgc_class)
            filenames = []
            for cutoff in cutoff_list:
                filenames.append(os.path.join(pathBase, "{}_c{:.2f}.network".format(bgc_class, cutoff)))
            cutoffs_and_filenames = list(zip(cutoff_list, filenames))
            del filenames[:]
            write_network_matrix(network_matrix, cutoffs_and_filenames, include_singletons, clusterNames, bgc_info)

            print("  Calling Gene Cluster Families")
            reduced_network = []
            pos_alignments = {}
            for row in network_matrix:
                reduced_network.append([int(row[0]), int(row[1]), row[2]])
                reverse = False
                if row[-1] == 1.0:
                    reverse = True
                pa = pos_alignments.setdefault(int(row[0]),{})
                # lcsStartA, lcsStartB, seedLength, reverse={True,False}
                pa[int(row[1])] = (int(row[-4]), int(row[-3]), int(row[-2]), reverse)
            del network_matrix[:]

            family_data = clusterJsonBatch(BGC_classes[bgc_class], pathBase, bgc_class,
                reduced_network, pos_alignments, clusterNames, bgc_info, 
                mibig_set, pfd_folder, bgc_fasta_folder, DomainList, 
                BGCs, AlignedDomainSequences, DomainCountGene, BGCGeneOrientation,
                cutoffs=cutoff_list, clusterClans=options.clans, clanCutoff=options.clan_cutoff, 
                htmlFolder=network_html_folder)
            for network_html_folder_cutoff in family_data:
                rundata_networks_per_run[network_html_folder_cutoff].append(family_data[network_html_folder_cutoff])
                if (len(family_data[network_html_folder_cutoff]["families"]) > 0):
                    html_subs_per_run[network_html_folder_cutoff].append({ "name" : bgc_class, "css" : bgc_class, "label" : bgc_class})
            del BGC_classes[bgc_class][:]
            del reduced_network[:]

    # fetch genome list for overview.js
    genomes = []
    classes = []
    clusterNamesToGenomes = {}
    clusterNamesToClasses = {}
    inputClustersIdx = [] # contain only indexes (from clusterNames) of input BGCs (non-mibig)
    for idx, bgc in enumerate(clusterNames):
        if bgc in mibig_set:
            continue
        inputClustersIdx.append(idx)
        # get class info
        product = bgc_info[bgc].product
        predicted_class = sort_bgc(product)
        if predicted_class not in classes:
            clusterNamesToClasses[bgc] = len(classes)
            classes.append(predicted_class)
        else:
            clusterNamesToClasses[bgc] = classes.index(predicted_class)
        # get identifier info
        identifier = ""
        if len(bgc_info[bgc].organism) > 1:
            identifier = bgc_info[bgc].organism
        else : # use original genome file name (i.e. exclude "..clusterXXX from antiSMASH run")
            file_name_base = os.path.splitext(os.path.basename(genbankDict[bgc][0]))[0]
            identifier = file_name_base.rsplit(".cluster",1)[0].rsplit(".region", 1)[0]
        if len(identifier) < 1:
            identifier = "Unknown Genome {}".format(len(genomes))
        if identifier not in genomes:
            clusterNamesToGenomes[bgc] = len(genomes)
            genomes.append(identifier)
        else:
            clusterNamesToGenomes[bgc] = genomes.index(identifier)
    run_data["input"]["accession"] = [{ "id": "genome_{}".format(i), "label": acc } for i, acc in enumerate(genomes)]
    run_data["input"]["accession_newick"] = [] # todo ...
    run_data["input"]["classes"] = [{ "label": cl } for cl in classes ] # todo : colors
    run_data["input"]["bgc"] = [{ "id": clusterNames[idx], "acc": clusterNamesToGenomes[clusterNames[idx]], "class": clusterNamesToClasses[clusterNames[idx]] } for idx in inputClustersIdx]


    # update family data (convert global bgc indexes into input-only indexes)
    for network_key in rundata_networks_per_run:
        for network in rundata_networks_per_run[network_key]:
            for family in network["families"]:
                new_members = []
                mibig = []
                for bgcIdx in family["members"]:
                    if bgcIdx in inputClustersIdx:                    
                        new_members.append(inputClustersIdx.index(bgcIdx))
                    else: # is a mibig bgc
                        clusterName = clusterNames[bgcIdx]
                        if clusterName in mibig_set:
                            mibig.append(clusterName)
                family["mibig"] = mibig
                family["members"] = new_members


    # generate overview data
    end_time = time.localtime()
    duration = int(time.mktime(end_time)) - int(time.mktime(start_time))
    run_data["end_time"] = time.strftime("%d/%m/%Y %H:%M:%S", end_time)
    run_data["duration"] = "{}h{}m{}s".format((duration // 3600), ((duration % 3600) // 60), ((duration % 3600) % 60))

    for cutoff in cutoff_list:
        # update overview.html
        html_folder_for_this_cutoff = "{}_c{:.2f}".format(network_html_folder, cutoff)
        run_data_for_this_cutoff = run_data.copy()
        run_data_for_this_cutoff["networks"] = rundata_networks_per_run[html_folder_for_this_cutoff]
        with open(os.path.join(html_folder_for_this_cutoff, "run_data.js"), "w") as run_data_js:
            run_data_js.write("var run_data={};\n".format(json.dumps(run_data_for_this_cutoff, indent=4, separators=(',', ':'), sort_keys=True)))
            run_data_js.write("dataLoaded();\n");
        # update bgc_results.js
        add_to_bigscape_results_js("{}_c{:.2f}".format(run_name, cutoff), html_subs_per_run[html_folder_for_this_cutoff], os.path.join(output_folder, "html_content", "js", "bigscape_results.js"))

    pickle.dump(bgc_info,open(os.path.join(cache_folder,'bgc_info.dict'),'wb'))
    runtime = time.time()-time1
    runtime_string = "\n\n\tMain function took {:.3f} s".format(runtime)
    with open(os.path.join(log_folder, "runtimes.txt"), 'a') as timings_file:
        timings_file.write(runtime_string + "\n")
    print(runtime_string)
    
