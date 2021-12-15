import os

from argparse import ArgumentParser
from multiprocessing import cpu_count


def CMD_parser():
    parser = ArgumentParser(prog="BiG-SCAPE")
    
    parser.add_argument("-l", "--label", dest="label", help="An extra label for\
                        this run (will be used as part of the folder name within\
                        the network_files results)")
    
    parser.add_argument("-i", "--inputdir", dest="inputdir", 
                        default=os.path.dirname(os.path.realpath(__file__)),
                        help="Input directory of gbk files, if left empty, all \
                        gbk files in current and lower directories will be used.")

    parser.add_argument("-o", "--outputdir", dest="outputdir", default="", 
                        required=True, help="Output directory, this will contain \
                        all output data files.")
    
    parser.add_argument("--pfam_dir", dest="pfam_dir",
                      default=os.path.dirname(os.path.realpath(__file__)), 
                      help="Location of hmmpress-processed Pfam files. Default \
                      is same location of BiG-SCAPE")
    
    parser.add_argument("-c", "--cores", dest="cores", default=cpu_count(),
                      help="Set the number of cores the script may use (default:\
                      use all available cores)")

    parser.add_argument("--include_gbk_str", dest="include_gbk_str",
                        default=['cluster', 'region'], nargs="+",
                        help="Only gbk files with this string(s) will be used \
                        for the analysis (default: 'cluster', 'region'). Use an \
                        asterisk to accept every file (overrides \
                        '--exclude_gbk_str')")

    parser.add_argument("--exclude_gbk_str", dest="exclude_gbk_str", 
                        default=['final'], nargs="+",
                        help="If any string in this list occurs in the gbk \
                        filename, this file will not be used for the analysis\
                         (default: final).")
        
    parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", 
                        default=False, help="Prints more detailed information. \
                        Toggle to activate.")
    
    parser.add_argument("--include_singletons", dest="include_singletons", 
                        action="store_true", default=False, help="Include nodes \
                        that have no edges to other nodes from the network. \
                        Toggle to activate.")
    
    parser.add_argument("-d", "--domain_overlap_cutoff", 
                        dest="domain_overlap_cutoff", default=0.1, help="Specify\
                        at which overlap percentage domains are considered to \
                        overlap. Domain with the best score is kept \
                        (default=0.1).")
    
    parser.add_argument("-m", "--min_bgc_size", dest="min_bgc_size", default=0,
                      help="Provide the minimum size of a BGC to be included in\
                      the analysis. Default is 0 base pairs")
    
    parser.add_argument("--mix", dest="mix", action="store_true", default=False, 
                        help="By default, BiG-SCAPE separates the analysis \
                        according to the BGC product (PKS Type I, NRPS, RiPPs, etc.) \
                        and will create network directories for each class. \
                        Toggle to include an analysis mixing all classes")
    
    parser.add_argument("--no_classify", dest="no_classify", action="store_true", 
                        default=False, help="By default, BiG-SCAPE classifies \
                        the output files analysis based on the BGC product. \
                        Toggle to deactivate (note that if the --mix parameter \
                        is not activated, BiG-SCAPE will not create any network \
                        file).")
    
    parser.add_argument("--banned_classes", nargs='+', dest="banned_classes", 
                        default=[], choices=["PKSI", "PKSother", "NRPS", "RiPPs", 
                                             "Saccharides", "Terpene", 
                                             "PKS-NRP_Hybrids", "Others"], 
                        help="Classes that should NOT be included in the \
                        classification. E.g. \"--banned_classes PKSI PKSOther\"")

    parser.add_argument("--cutoffs", dest="cutoffs", nargs="+", default=[0.30], 
                        type=float, help="Generate networks using multiple raw \
                        distance cutoff values. Values should be in the range \
                        [0.0, 1.0]. Example: --cutoffs 0.1 0.25 0.5 1.0. Default: \
                        c=0.3.")
                    
    parser.add_argument("--clans-off", dest="clans", action="store_false", 
                        default=True, help="Toggle to deactivate a second \
                        layer of clustering to attempt to group families \
                        into clans")

    parser.add_argument("--clan_cutoff",dest="clan_cutoff",default=[0.3,0.7], 
                        type=float, nargs=2, help="Cutoff Parameters for which \
                        clustering families into clans will be performed in raw \
                        distance. First value is the cutoff value family \
                        assignments for BGCs used in clan clustering (default: \
                        0.3). Second value is the cutoff value for clustering \
                        families into clans (default: 0.7). Average linkage for \
                        BGCs in a family is used for distances between families. \
                        Valid values are in the range [0.0, 1.0]. Example: \
                        --clan_cutoff 0.3 0.7)")

    parser.add_argument("--hybrids-off", dest="hybrids", action="store_false", 
                        default=True, help="Toggle to also add BGCs with hybrid\
                        predicted products from the PKS/NRPS Hybrids and Others\
                        classes to each subclass (e.g. a 'terpene-nrps' BGC from\
                        Others would be added to the Terpene and NRPS classes)")
    
    parser.add_argument("--mode", dest="mode", default="glocal", choices=["global",
                            "glocal", "auto"], help="Alignment mode for each pair of\
                            gene clusters. 'global': the whole list of domains \
                            of each BGC are compared; 'glocal': Longest Common \
                            Subcluster mode. Redefine the subset of the domains \
                            used to calculate distance by trying to find the \
                            longest slice of common domain content per gene in \
                            both BGCs, then expand each slice. 'auto': use glocal\
                            when at least one of the BGCs in each pair has the \
                            'contig_edge' annotation from antiSMASH v4+, \
                            otherwise use global mode on that pair")
    
    parser.add_argument("--anchorfile", dest="anchorfile", 
                        default=os.path.join(os.path.dirname(os.path.realpath(__file__)),"anchor_domains.txt"),
                        help="Provide a custom location for the anchor domains \
                        file, default is anchor_domains.txt.")
                    
    parser.add_argument("--force_hmmscan", dest="force_hmmscan", action="store_true", 
                        default=False, help="Force domain prediction using \
                        hmmscan even if BiG-SCAPE finds processed domtable files\
                        (e.g. to use a new version of PFAM).")
    
    parser.add_argument("--skip_ma", dest="skip_ma", action="store_true", 
                        default=False, help="Skip multiple alignment of domains'\
                        sequences. Use if alignments have been generated in a \
                        previous run.")
    
    parser.add_argument("--mibig", dest="mibig21", default=False, action="store_true",
                        help="Use included BGCs from then MIBiG database. Only \
                        relevant (i.e. those with distance < max(cutoffs) against\
                        the input set) will be used. Currently uses version 2.1\
                        of MIBiG. See https://mibig.secondarymetabolites.org/")
    
    parser.add_argument("--mibig14", dest="mibig14", default=False, action="store_true",
                        help="Include BGCs from version 1.4 of MIBiG")
    
    parser.add_argument("--mibig13", dest="mibig13", default=False, action="store_true",
                        help="Include BGCs from version 1.3 of MIBiG")
    
    parser.add_argument("--query_bgc", help="Instead of making an all-VS-all \
                        comparison of all the input BGCs, choose one BGC to \
                        compare with the rest of the set (one-VS-all). The \
                        query BGC does not have to be within inputdir")
    
    parser.add_argument("--domain_includelist", help="Only analyze BGCs that \
                        include domains with the pfam accessions found in the \
                        domain_includelist.txt file", default=False,
                        action="store_true")

    parser.add_argument("--version", action="version", version="%(prog)s 1.1.2 (2021-06-03)")

    return parser.parse_args()
