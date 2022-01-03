import os
import sys

from Bio import SeqIO
from itertools import combinations

from bgctools.bgctools import sort_bgc

def process_gbk_files(
        gbk,
        min_bgc_size,
        bgc_info,
        files_no_proteins,
        files_no_biosynthetic_genes,
        verbose,
        bgc_fasta_folder,
        force_hmmscan,
        valid_classes,
        bgc_data,
        genbankDict,

):
    """ Given a file path to a GenBank file, reads information about the BGC"""

    biosynthetic_genes = set()
    product_list_per_record = []
    fasta_data = []
    save_fasta = True
    adding_sequence = False
    contig_edge = False
    total_seq_length = 0
    record_end = 0
    offset_record_position = 0
    bgc_locus_tags = []
    locus_sequences = {}
    locus_coordinates = {}

    file_folder, fname = os.path.split(gbk)
    clusterName = fname[:-4]

    # See if we need to keep the sequence
    # (Currently) we have to open the file anyway to read all its 
    # properties for bgc_info anyway...
    outputfile = os.path.join(bgc_fasta_folder, clusterName + '.fasta')
    if os.path.isfile(outputfile) and os.path.getsize(outputfile) > 0 and not force_hmmscan:
        if verbose:
            print(" File {} already processed".format(outputfile))
        save_fasta = False
    else:
        save_fasta = True
    
    try:
        # basic file verification. Substitutes check_data_integrity
        records = list(SeqIO.parse(gbk, "genbank"))
    except ValueError as e:
        print("   Error with file {}: \n    '{}'".format(gbk, str(e)))
        print("    (This file will be excluded from the analysis)")
        return
    else:
        total_seq_length = 0
        bgc_size = 0
        cds_ctr = 0
        product = "no type"
        offset_record_position = 0
        
        max_width = 0 # This will be used for the SVG figure
        record_count = 0
        
        for record in records:
            record_count += 1
            bgc_size += len(record.seq)
            if len(record.seq) > max_width:
                max_width = len(record.seq)
            
            for feature in record.features:
                # antiSMASH <= 4
                if feature.type == "cluster":
                    if "product" in feature.qualifiers:
                        # in antiSMASH 4 there should only be 1 product qualifiers
                        for product in feature.qualifiers["product"]:
                            for p in product.replace(" ","").split("-"):
                                product_list_per_record.append(p)
                                
                    if "contig_edge" in feature.qualifiers:
                        # there might be mixed contig_edge annotations
                        # in multi-record files. Turn on contig_edge when
                        # there's at least one annotation
                        if feature.qualifiers["contig_edge"][0] == "True":
                            if verbose:
                                print(" Contig edge detected in {}".format(fname))
                            contig_edge = True
                        
                # antiSMASH = 5
                if "region" in feature.type:
                    if "product" in feature.qualifiers:
                        for product in feature.qualifiers["product"]:
                            product_list_per_record.append(product)
                            
                    if "contig_edge" in feature.qualifiers:
                        # there might be mixed contig_edge annotations
                        # in multi-record files. Turn on contig_edge when
                        # there's at least one annotation
                        if feature.qualifiers["contig_edge"][0] == "True":
                            if verbose:
                                print(" Contig edge detected in {}".format(fname))
                            contig_edge = True
                            
                        
                # Get biosynthetic genes + sequences
                if feature.type == "CDS":
                    cds_ctr += 1
                    CDS = feature
                    
                    gene_id = ""
                    if "gene" in CDS.qualifiers:
                        gene_id = CDS.qualifiers.get('gene',"")[0]
                        
                    
                    protein_id = ""
                    if "protein_id" in CDS.qualifiers:
                        protein_id = CDS.qualifiers.get('protein_id',"")[0]
                    
                    # nofuzzy_start/nofuzzy_end are obsolete
                    # http://biopython.org/DIST/docs/api/Bio.SeqFeature.FeatureLocation-class.html#nofuzzy_start
                    gene_start = offset_record_position + max(0, int(CDS.location.start))
                    gene_end = offset_record_position + max(0, int(CDS.location.end))
                    record_end = gene_end
                    
                    direction = CDS.location.strand
                    if direction == 1:
                        strand = '+'
                    else:
                        strand = '-'
                        
                    fasta_header = "{}_ORF{}:gid:{}:pid:{}:loc:{}:{}:strand:{}".format(clusterName, str(cds_ctr), str(gene_id).replace(":","_"), str(protein_id).replace(":","_"), str(gene_start), str(gene_end), strand)
                    fasta_header = fasta_header.replace(">","") #the coordinates might contain larger than signs, tools upstream don't like this
                    fasta_header = fasta_header.replace(" ", "") #the domtable output format (hmmscan) uses spaces as a delimiter, so these cannot be present in the fasta header

                    # antiSMASH <=4
                    if "sec_met" in feature.qualifiers:
                        if "Kind: biosynthetic" in feature.qualifiers["sec_met"]:
                            biosynthetic_genes.add(fasta_header)

                    # antiSMASH == 5
                    if "gene_kind" in feature.qualifiers:
                        if "biosynthetic" in feature.qualifiers["gene_kind"]:
                            biosynthetic_genes.add(fasta_header)
                    
                    fasta_header = ">"+fasta_header
                    

                    if 'translation' in CDS.qualifiers.keys():
                        prot_seq = CDS.qualifiers['translation'][0]
                    # If translation isn't available translate manually, this will take longer
                    else:
                        nt_seq = CDS.location.extract(record.seq)
                        
                        # If we know sequence is an ORF (like all CDSs), codon table can be
                        #  used to correctly translate alternative start codons.
                        #  see http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc25
                        # If the sequence has a fuzzy start/end, it might not be complete,
                        # (therefore it might not be the true start codon)
                        # However, in this case, if 'translation' not available, assume 
                        #  this is just a random sequence 
                        complete_cds = False 
                        
                        # More about fuzzy positions
                        # http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc39
                        fuzzy_start = False 
                        if str(CDS.location.start)[0] in "<>":
                            complete_cds = False
                            fuzzy_start = True
                            
                        fuzzy_end = False
                        if str(CDS.location.end)[0] in "<>":
                            fuzzy_end = True
                        
                        #for protein sequence if it is at the start of the entry assume 
                        # that end of sequence is in frame and trim from the beginning
                        #if it is at the end of the genbank entry assume that the start 
                        # of the sequence is in frame
                        reminder = len(nt_seq)%3
                        if reminder > 0:
                            if fuzzy_start and fuzzy_end:
                                print("Warning, CDS ({}, {}) has fuzzy\
                                    start and end positions, and a \
                                    sequence length not multiple of \
                                    three. Skipping".format(clusterName, 
                                    CDS.qualifiers.get('locus_tag',"")[0]))
                                break
                            
                            if fuzzy_start:
                                if reminder == 1:
                                    nt_seq = nt_seq[1:]
                                else:
                                    nt_seq = nt_seq[2:]
                            # fuzzy end
                            else:
                                #same logic reverse direction
                                if reminder == 1:
                                    nt_seq = nt_seq[:-1]
                                else:
                                    nt_seq = nt_seq[:-2]
                        
                        # The Genetic Codes: www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
                        if "transl_table" in CDS.qualifiers.keys():
                            CDStable = CDS.qualifiers.get("transl_table", "")[0]
                            prot_seq = str(nt_seq.translate(table=CDStable, to_stop=True, cds=complete_cds))
                        else:
                            prot_seq = str(nt_seq.translate(to_stop=True, cds=complete_cds))
                            
                    total_seq_length += len(prot_seq)
                
                
                    bgc_locus_tags.append(fasta_header)
                    locus_sequences[fasta_header] = prot_seq
                    locus_coordinates[fasta_header] = (gene_start, gene_end, len(prot_seq))
                    

            # TODO: if len(biosynthetic_genes) == 0, traverse record again
            # and add CDS with genes that contain domains labeled sec_met
            # we'll probably have to have a list of domains if we allow
            # fasta files as input
            
            # make absolute positions for ORFs in next records
            offset_record_position += record_end + 100
        
        if bgc_size > min_bgc_size:  # exclude the bgc if it's too small
            # check what we have product-wise
            # In particular, handle different products for multi-record files
            product_set = set(product_list_per_record)
            if len(product_set) == 1: # only one type of product
                product = product_list_per_record[0]
            elif "other" in product_set: # more than one, and it contains "other"
                if len(product_set) == 2:
                    product = list(product_set - {'other'})[0] # product = not "other"
                else:
                    product = ".".join(product_set - {'other'}) # likely a hybrid
            else:
                product = ".".join(product_set) # likely a hybrid
                
            # Don't keep this bgc if its type not in valid classes specified by user
            # This will avoid redundant tasks like domain detection
            subproduct = set()
            for p in product.split("."):
                subproduct.add(sort_bgc(p).lower())
            if "nrps" in subproduct and ("pksi" in subproduct or "pksother" in subproduct):
                subproduct.add("pks-nrp_hybrids")
                
            if len(valid_classes & subproduct) == 0:
                if verbose:
                    print(" Skipping {} (type: {})".format(clusterName, product))
                return False
            
            # assuming that the definition field is the same in all records
            # product: antiSMASH predicted class of metabolite
            # gbk definition
            # number of records (for Arrower figures)
            # max_width: width of the largest record (for Arrower figures)
            # id: the GenBank's accession
            #bgc_info[clusterName] = (product, records[0].description, len(records), max_width, records[0].id, biosynthetic_genes.copy())
            # TODO contig_edge annotation is not present for antiSMASH v < 4
            # Perhaps we can try to infer if it's in a contig edge: if
            # - first biosynthetic gene start < 10kb or
            # - max_width - last biosynthetic gene end < 10kb (but this will work only for the largest record)
            bgc_info[clusterName] = bgc_data(records[0].id, records[0].description, product, len(records), max_width, bgc_size + (record_count-1)*1000, records[0].annotations["organism"], ",".join(records[0].annotations["taxonomy"]), biosynthetic_genes.copy(), contig_edge)

            if len(bgc_info[clusterName].biosynthetic_genes) == 0:
                files_no_biosynthetic_genes.append(clusterName+".gbk")

            # TODO why re-process everything if it was already in the list?
            # if name already in genbankDict.keys -> add file_folder
            # else: extract all info
            if clusterName in genbankDict.keys():
                # Name was already in use. Use file_folder as the new sample's name
                genbankDict[clusterName][1].add(file_folder) 
            else:
                # See if we need to write down the sequence
                if total_seq_length > 0:
                    # location of first instance of the file is genbankDict[clustername][0]
                    genbankDict.setdefault(clusterName, [gbk, set([file_folder])])

                    if save_fasta:
                        # Find overlaps in CDS regions and delete the shortest ones.
                        # This is thought as a solution for selecting genes with 
                        # alternate splicing events
                        # Food for thought: imagine CDS A overlapping CDS B overlapping
                        # CDS C. If len(A) > len(B) > len(C) and we first compare A vs B
                        # and delete A, then B vs C and delete B: would that be a better
                        # solution than removing B? Could this actually happen?
                        # TODO What if the overlapping CDS is in the reverse strand?
                        #  maybe it should be kept as it is
                        # TODO what are the characterized differences in prokarytote
                        #  vs eukaryote CDS overlap?
                        del_list = set()
                        for a, b in combinations(bgc_locus_tags, 2):
                            a_start, a_end, a_len = locus_coordinates[a]
                            b_start, b_end, b_len = locus_coordinates[b]
                            
                            if b_end <= a_start or b_start >= a_end:
                                pass
                            else:
                                # calculate overlap
                                if a_start > b_start:
                                    ov_start = a_start
                                else:
                                    ov_start = b_start

                                if a_end < b_end:
                                    ov_end = a_end
                                else:
                                    ov_end = b_end

                                overlap_length = ov_end - ov_start
                                
                                # allow the overlap to be as large as 10% of the
                                # shortest CDS. Overlap length is in nucleotides
                                # here, whereas a_len, b_len are protein 
                                # sequence lengths
                                if overlap_length/3 > 0.1*min(a_len,b_len):
                                    if a_len > b_len:
                                        del_list.add(b)
                                    else:
                                        del_list.add(a)
                        
                        for locus in del_list:
                            if verbose:
                                print("   Removing {} because it overlaps with other ORF".format(locus))
                            bgc_locus_tags.remove(locus)
                        
                        with open(outputfile,'w') as fastaHandle:
                            for locus in bgc_locus_tags:
                                fastaHandle.write("{}\n".format(locus))
                                fastaHandle.write("{}\n".format(locus_sequences[locus]))
                            adding_sequence = True
                else:
                    files_no_proteins.append(fname)

            if verbose:
                print("  Adding {} ({} bps)".format(fname, str(bgc_size)))
                                
        else:
            print(" Discarding {} (size less than {} bp, was {})".format(clusterName, str(min_bgc_size), str(bgc_size)))
    
    return adding_sequence


def get_gbk_files(
        inputpath, outputdir, bgc_fasta_folder, min_bgc_size, include_gbk_str,
        exclude_gbk_str, bgc_info, mode, verbose, force_hmmscan, valid_classes, bgc_data, genbankDict
):
    """Searches given directory for genbank files recursively, will assume that
    the genbank files that have the same name are the same genbank file. 
    Returns a dictionary that contains the names of the clusters found as keys
    and a list that contains [0] a path to the genbank file and [1] the 
    samples that the genbank file is a part of.
    Extract and write the sequences as fasta files if not already in the Fasta 
    folder.
    return: {cluster_name:[genbank_path,[s_a,s_b...]]}
    """
    file_counter = 0
    processed_sequences = 0
    files_no_proteins = []
    files_no_biosynthetic_genes = []


    if os.path.isfile(inputpath):
        files = [inputpath]
    else:
        # Unfortunately, this does not work in Python 2:
        #files = glob(os.path.join(inputpath,"**/*.gbk"), recursive=True) 
        files = [os.path.join(dirpath, f) for dirpath, dirnames, files in os.walk(inputpath)
                 for f in files if f.endswith(".gbk")]
        
    for filepath in files:
        file_folder, fname = os.path.split(filepath)
        
        if len(include_gbk_str) == 1 and include_gbk_str[0] == "*":
            pass
        else:
            if not any([word in fname for word in include_gbk_str]):
                continue
            
            if exclude_gbk_str != [] and any([word in fname for word in exclude_gbk_str]):
                print(" Skipping file " + fname)
                continue
        
        if "_ORF" in fname:
            print(" Skipping file {} (string '_ORF' is used internally)".format(fname))
            continue
        
        if " " in filepath:
            sys.exit("\nError: Input GenBank files should not have spaces in their path as hmmscan cannot process them properly ('too many arguments').")
        
        file_counter += 1
        if process_gbk_files(filepath, min_bgc_size, bgc_info, files_no_proteins, files_no_biosynthetic_genes, verbose, bgc_fasta_folder, force_hmmscan, valid_classes, bgc_data, genbankDict):
            processed_sequences += 1
    
    if len(files_no_proteins) > 0:
        print("  Warning: Input set has files without protein sequences. They will be discarded")
        print("   (See no_sequences_list.txt)")
        with open(os.path.join(outputdir, "no_sequences_list.txt"), "w") as noseqs:
            for f in sorted(files_no_proteins):
                noseqs.write("{}\n".format(f))
        
    if len(files_no_biosynthetic_genes) > 0 and (mode == "glocal" or mode == "auto"):
        print("  Warning: Input set has files with no Biosynthetic Genes (affects alignment mode)")
        print("   See no_biosynthetic_genes_list.txt")
        with open(os.path.join(outputdir, "logs", "no_biosynthetic_genes_list.txt"), "w") as nobiogenes:
            for f in sorted(files_no_biosynthetic_genes):
                nobiogenes.write("{}\n".format(f))
    
    print("\n Starting with {:d} files".format(file_counter))
    print(" Files that had its sequence extracted: {:d}".format(processed_sequences))

    return
