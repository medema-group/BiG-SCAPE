import logging
import os
import sys
import numpy as np
import networkx as nx
import subprocess
import json
import shutil


from collections import defaultdict
from operator import itemgetter
from itertools import combinations
from sklearn.cluster import AffinityPropagation
from scipy.optimize import linear_sum_assignment
from Bio import Phylo

from src.bgctools import get_composite_bgc_similarities
from src.big_scape.bgc_collection import BgcCollection
from src.utility.io import create_directory


def clusterJsonBatch(bgcs, pathBase, className, matrix, pos_alignments, bgc_collection: BgcCollection,
                     mibig_set, pfd_folder, bgc_fasta_folder,
                     AlignedDomainSequences, cutoffs=[1.0],
                     damping=0.9, clusterClans=False, clanCutoff=(0.5, 0.8), htmlFolder=None,
                     verbose=False):
    """BGC Family calling
    Uses csr sparse matrices to call Gene Cluster Families (GCFs) using Affinity
    Propagation.

    Cutoff values are in distance (i.e. if two clusters are further than cutoff
    value, similarity is 0)
    Larger cutoff values are more permissive

    bgcs: ordered list of integers (ascending, unique, but not necessarily
        consecutive) representing the index in the main clusterNames list. Every
        number here is an "external" index
    matrix: list of lists of idx0, idx1, d where the first two elements correspond
        to indices from `bgcs`
    pathBase: folder where GCF files will be deposited
    """
    numBGCs = len(bgcs)

    simDict = {} # dictionary of dictionaries
    # Doing this so it only has to go through the matrix once
    for row in matrix:
        gc1, gc2, distance = row

        if distance < 1.0:
            similarity = 1 - distance
        else:
            similarity = 0
        gcSimilarities = simDict.setdefault(gc1, {})
        gcSimilarities[gc2] = similarity

    clanClassificationCutoff, clanDistanceCutoff = clanCutoff
    if clusterClans:
        logging.debug('Clustering Clans Enabled with parameters clanClassificationCutoff: %s, \
                      clanDistanceCutoff: %s', clanClassificationCutoff, clanDistanceCutoff)

    # External index number to Internal, consecutive, index
    # If needed, bgcInt2Ext[i] would simply be bgcs[i]
    bgcExt2Int = dict(zip(bgcs, range(numBGCs)))

    # Get info on all BGCs to export to .js for visualization
    bs_data = []
    bgcJsonDict = {}
    for bgc in bgcs:
        bgcName = bgc_collection.bgc_name_tuple[bgc]
        bgcJsonDict[bgcName] = {}
        bgcJsonDict[bgcName]["id"] = bgcName
        bgcJsonDict[bgcName]["desc"] = bgc_collection.bgc_collection_dict[bgcName].bgc_info.description
        bgcJsonDict[bgcName]["start"] = 1
        bgcJsonDict[bgcName]["end"] = bgc_collection.bgc_collection_dict[bgcName].bgc_info.bgc_size
        bgcJsonDict[bgcName]["mibig"] = bgcName in mibig_set

        pfdFile = os.path.join(pfd_folder, bgcName + ".pfd")
        fastaFile = os.path.join(bgc_fasta_folder, bgcName + ".fasta")

        orfDict = defaultdict(dict)

        ## read fasta file first to get orfs
        # We cannot get all the info exclusively from the pfd because that only
        # contains ORFs with predicted domains (and we need to draw empty genes
        # as well)
        for line in open(fastaFile):
            if line[0] == ">":
                header = line.strip()[1:].split(':')
                orf = header[0]
                if header[2]:
                    orfDict[orf]["id"] = header[2]
                elif header[4]:
                    orfDict[orf]["id"] = header[4]
                else:
                    orfDict[orf]["id"] = orf

                ## broken gene goes into cluster, need this so js doesn't throw an error
                if int(header[6]) <= 1:
                    orfDict[orf]["start"] = 1
                else:
                    orfDict[orf]["start"] = int(header[6])

                orfDict[orf]["end"] = int(header[7])

                if header[-1] == '+':
                    orfDict[orf]["strand"] = 1
                else:
                    orfDict[orf]["strand"] = -1

                orfDict[orf]["domains"] = []

        ## now read pfd file to add the domains to each of the orfs
        for line in open(pfdFile):
            entry = line.split('\t')
            orf = entry[-1].strip().split(':')[0]
            pfamID = entry[5].split('.')[0]
            orfDict[orf]["domains"].append({'code': pfamID, 'start': int(entry[3]), 'end': int(entry[4]), 'bitscore': float(entry[1])})
        # order of ORFs is important here because I use it to get a translation
        # between the "list of ORFs with domains" to "list of all ORFs" later on
        bgcJsonDict[bgcName]['orfs'] = sorted(orfDict.values(), key=itemgetter("start"))
    bs_data = [bgcJsonDict[bgc_collection.bgc_name_tuple[bgc]] for bgc in bgcs]


    # Create network
    g = nx.Graph()

    results = {}
    for cutoff in cutoffs:
        family_data = { # will be returned, for use in overview.js
            "label": className,
            "families": [],
            "families_similarity": []
        }

        logging.info("  Cutoff: %s", cutoff)
        # first task is to find labels (exemplars) for each node.
        # Assign disconnected (singleton) BGCs to themselves
        labels = [0 for i in range(numBGCs)]

        # clear all edges from network
        g.clear()

        # add all nodes
        g.add_nodes_from(bgcs)

        # add all (allowed) edges
        for bgc1 in bgcs:
            for bgc2 in simDict.get(bgc1,{}).keys():
                if simDict[bgc1][bgc2] > 1 - cutoff:
                    g.add_edge(bgc1, bgc2)

        for subgraph in nx.connected_components(g):
            numBGCs_subgraph = len(subgraph)
            # smaller subgraphs don't need to be clustered
            if numBGCs_subgraph < 3:
                temp = list(subgraph)
                for bgc in temp:
                    labels[bgcExt2Int[bgc]] = bgcExt2Int[temp[0]]
            else:
                bgcExt2Sub_ = dict(zip([c for c in subgraph], range(numBGCs_subgraph)))
                bgcSub2Ext_ = dict(zip(range(numBGCs_subgraph), [c for c in subgraph]))

                simMatrix = np.zeros((numBGCs_subgraph, numBGCs_subgraph), dtype=np.float32)
                for bgc1, bgc2 in combinations(subgraph,2):
                    try:
                        simMatrix[bgcExt2Sub_[bgc1], bgcExt2Sub_[bgc2]] = simDict[bgc1][bgc2]
                    except KeyError:
                        simMatrix[bgcExt2Sub_[bgc1], bgcExt2Sub_[bgc2]] = simDict[bgc2][bgc1]

                af = AffinityPropagation(damping=damping, max_iter=1000, convergence_iter=200, affinity="precomputed").fit(simMatrix)
                labelsSub = af.labels_
                exemplarsSub = af.cluster_centers_indices_

                # TODO go straight to familiesDict
                for i in range(numBGCs_subgraph):
                    labels[bgcExt2Int[bgcSub2Ext_[i]]] = bgcExt2Int[bgcSub2Ext_[exemplarsSub[labelsSub[i]]]]

        # Recalculate distance matrix as we'll need it with clans
        #simMatrix = lil_matrix((numBGCs, numBGCs), dtype=np.float32)
        simMatrix = np.zeros((numBGCs, numBGCs), dtype=np.float32)
        for bgc1 in bgcs:
            # first make sure it is similar to itself
            simMatrix[bgcExt2Int[bgc1],bgcExt2Int[bgc1]] = 1
            for bgc2 in simDict.get(bgc1,{}).keys():
                # you might get 0 values if there were matrix entries under the
                # cutoff. No need to input these in the sparse matrix

                if simDict[bgc1][bgc2] > 1-cutoff:
                    # Ensure symmetry
                    simMatrix[bgcExt2Int[bgc1], bgcExt2Int[bgc2]] = simDict[bgc1][bgc2]
                    simMatrix[bgcExt2Int[bgc2], bgcExt2Int[bgc1]] = simDict[bgc1][bgc2]

        logging.debug("   ...done")

        bs_distances = [[float("{:.3f}".format(simMatrix[row, col])) for col in
                         range(row+1)] for row in range(numBGCs)]


        familiesDict = defaultdict(list)
        for i in range(numBGCs):
            familiesDict[bgcs[labels[i]]].append(bgcs[i])
        familyIdx = sorted(familiesDict.keys()) # identifiers for each family


        ##
        ## Get conserved domain core information to build phylogenetic tree
        ##
        gcf_trees_path = os.path.join(pathBase, "GCF_trees")
        if not os.path.exists(gcf_trees_path):
            os.makedirs(gcf_trees_path)

        newick_trees = {} # key:original bgc index
        for exemplar_idx in familiesDict:
            exemplar = bgc_collection.bgc_name_tuple[exemplar_idx]
            gcf = familiesDict[exemplar_idx][:]
            if len(gcf) < 3:
                newick_trees[exemplar_idx] = "({}):0.01;".format(",".join([str(bgcExt2Int[x])+":0.0" for x in gcf]))

                #print(newick_trees[exemplar_idx])
                #TODO make some default alignment data to send to the json file
                continue

            domain_sets = {}

            # make a frequency table (not counting copies):
            frequency_table = defaultdict(int)
            for bgc in gcf:
                domain_sets[bgc] = set(bgc_collection.bgc_ordered_domain_list[bgc_collection.bgc_name_tuple[bgc]])
                for domain in domain_sets[bgc]:
                    frequency_table[domain] += 1

            # Remove all PKS/NRPS domains
            # TODO but this was intended to be done when we were considering
            # directly using Corason. Should these domains still be removed?
            # We are now able to include them accurately with the Hungarian
            # matching
            # If implemented, fill nrps_pks_domains first!
            #domains_in_gcf = set(frequency_table.keys())
            #for erase_domain in (nrps_pks_domains & domains_in_gcf):
                #del frequency_table[erase_domain]

            # Find the set of [(tree domains)]. They should 1) be in the exemplar
            # and 2) appear with the most frequency. Iterate over the different
            # frequencies (descending) until set is not empty
            tree_domains = set()
            frequencies = sorted(set(frequency_table.values()), reverse=True)

            # first try with domain(s) with max frequency, even if it's just one
            f = 0
            while len(tree_domains) == 0 and f < len(frequencies):
                for domain in frequency_table:
                    if frequency_table[domain] == frequencies[f] and domain in domain_sets[exemplar_idx]:
                        tree_domains.add(domain)

                f += 1


            if len(tree_domains) == 1:
                logging.debug("  Note: core shared domains for GCF %s consists of a single domain (%s)", exemplar_idx, [x for x in tree_domains][0])

            # Get the alignments of the core domains
            alignments = {}
            # initialize every sequence alignment entry. Don't do defaultdict!
            alignments[exemplar_idx] = ""

            out_of_tree_bgcs = [] # bgcs that don't share a common domain core
            delete_list = set() # remove this bgcs from alignment
            gcf.remove(exemplar_idx) # separate exemplar from the rest of the bgcs
            for bgc in gcf:
                alignments[bgc] = ""

            match_dict = {}
            for domain in tree_domains:
                specific_domain_list_A = bgc_collection.bgc_collection_dict[exemplar].domain_name_info[domain]
                num_copies_a = len(specific_domain_list_A)
                for exemplar_domain_copy in specific_domain_list_A:
                    alignments[exemplar_idx] += AlignedDomainSequences[exemplar_domain_copy]

                seq_length = len(AlignedDomainSequences[specific_domain_list_A[0]])

                for bgc in alignments:
                    match_dict.clear()
                    if bgc == exemplar_idx:
                        pass
                    elif domain not in domain_sets[bgc]:
                        logging.debug("   BGC %s (%s) does not share a common domain core \
                                      (GCF: %s, domain: %s)", 
                                      bgc_collection.bgc_name_tuple[bgc], bgc, exemplar_idx, domain)
                        out_of_tree_bgcs.append(bgc)
                        delete_list.add(bgc)
                    elif bgc in delete_list:
                        pass
                    else:
                        specific_domain_list_B = bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].domain_name_info[domain]
                        num_copies_b = len(specific_domain_list_B)

                        DistanceMatrix = np.ndarray((num_copies_a,num_copies_b))

                        for domsa in range(num_copies_a):
                            for domsb in range(num_copies_b):
                                # TODO NOT taking into consideration any LCS slicing
                                # i.e. we're comparing ALL copies of this domain
                                sequence_tag_a = specific_domain_list_A[domsa]
                                sequence_tag_b = specific_domain_list_B[domsb]

                                aligned_seqA = AlignedDomainSequences[sequence_tag_a]
                                aligned_seqB = AlignedDomainSequences[sequence_tag_b]

                                matches = 0
                                gaps = 0

                                for position in range(seq_length):
                                    if aligned_seqA[position] == aligned_seqB[position]:
                                        if aligned_seqA[position] != "-":
                                            matches += 1
                                        else:
                                            gaps += 1

                                DistanceMatrix[domsa][domsb] = 1 - ( matches/(seq_length-gaps) )

                        BestIndexes = linear_sum_assignment(DistanceMatrix)
                        # at this point is not ensured that we have the same order
                        # for the exemplar's copies (rows in BestIndexes)
                        # ideally they should go from 0-numcopies. Better make sure
                        for x in range(len(BestIndexes[0])):
                            match_dict[BestIndexes[0][x]] = BestIndexes[1][x]

                        for copy in range(num_copies_a):
                            try:
                                alignments[bgc] += AlignedDomainSequences[specific_domain_list_B[match_dict[copy]]]
                            except KeyError:
                                # This means that this copy of exemplar did not
                                # have a match in bgc (i.e. bgc has less copies
                                # of this domain than exemplar)
                                alignments[bgc] += "-"*seq_length

                for bgc in list(delete_list):
                    del alignments[bgc]
                delete_list.clear()

            # need this to change the labels in the trees that are read from files
            bgc_name_to_idx = {}

            # save compiled alignments of the GCF domain core as fastas
            alignment_file_path = os.path.join(gcf_trees_path,"GCF_c{:4.2f}_{:05d}_alignment.fasta".format(cutoff,exemplar_idx))
            with open(alignment_file_path, "w") as gcf_alignment_file:
                gcf_alignment_file.write(">{}\n{}\n".format(exemplar, alignments[exemplar_idx]))
                bgc_name_to_idx[exemplar] = exemplar_idx

                for bgc in alignments:
                    if bgc != exemplar_idx:
                        gcf_alignment_file.write(">{}\n{}\n".format(bgc_collection.bgc_name_tuple[bgc], alignments[bgc]))
                        bgc_name_to_idx[bgc_collection.bgc_name_tuple[bgc]] = bgc

            # launch fasttree to make tree
            logging.debug("  Working GCF %s, cutoff %s", exemplar_idx, cutoff)

            # make tree
            newick_file_path = os.path.join(gcf_trees_path, "GCF_c{:4.2f}_{:05d}.newick".format(cutoff,exemplar_idx))
            with open(newick_file_path, "w") as newick_file:
                command = ["fasttree", "-nopr", "-quiet", alignment_file_path]
                p = subprocess.Popen(command, stdout=newick_file, shell=False)
                p.wait() # only with process has terminated will the file be ready

            # read tree, post-process it and save it
            if not os.path.isfile(newick_file_path) or os.path.getsize(newick_file_path) == 0:
                logging.error(newick_file_path)
                logging.error("The above newick file was not created or is empty (GCF_c%2f_%5d)", 
                              cutoff, exemplar_idx)
                sys.exit()
            else:
                with open(newick_file_path,"r") as newick_file:
                    try:
                        tree = Phylo.read(newick_file, 'newick')
                    except ValueError as e:
                        logging.warning(" There was an error while reading tree file %s", newick_file)
                        logging.warning(str(e))
                        newick_trees[exemplar_idx] = ""
                    else:
                        try:
                            tree.root_at_midpoint()
                        except UnboundLocalError:
                            # Noticed this could happen if the sequences are exactly
                            # the same and all distances == 0
                            
                            logging.debug(" Note: Unable to root at midpoint file %s", newick_file_path)
                            pass
                        newick = tree.format("newick")

                        # convert branches' names to indices for visualization
                        for name in bgc_name_to_idx:
                            newick = newick.replace(name, str(bgcExt2Int[bgc_name_to_idx[name]]))

                        newick_trees[exemplar_idx] = newick

        ### - - - GCC - - -
        bs_similarity_families = []
        if clusterClans and cutoff == clanClassificationCutoff:
            # Detect if there's only 1 GCF. It makes pySAPC crash
            if len(familyIdx) == 1:
                clanLabels = [1]
            else:
                #famSimMatrix = lil_matrix((len(familyIdx), len(familyIdx)), dtype=np.float32)
                famSimMatrix = np.zeros((len(familyIdx), len(familyIdx)), dtype=np.float32)
                familiesExt2Int = {gcfExtIdx:gcfIntIdx for gcfIntIdx,gcfExtIdx in enumerate(familyIdx)}

                for familyI, familyJ in [tuple(sorted(combo)) for combo in combinations(familyIdx, 2)]:
                    famSimilarities = []
                    # currently uses the average distance of all average distances
                    # between bgc from gcf I to all bgcs from gcf J
                    for bgcI in familiesDict[familyI]:
                        similarities = [simMatrix[bgcExt2Int[bgcI], bgcExt2Int[bgcJ]] for bgcJ in familiesDict[familyJ]]
                        famSimilarities.append(sum(similarities, 0.0) / len(similarities))

                    try:
                        familySimilarityIJ = sum(famSimilarities, 0.0)/len(famSimilarities)
                    except ZeroDivisionError:
                        familySimilarityIJ = 0.0

                    if familySimilarityIJ > 1 - clanDistanceCutoff:
                        # Ensure symmetry
                        famSimMatrix[familiesExt2Int[familyI], familiesExt2Int[familyJ]] = familySimilarityIJ
                        famSimMatrix[familiesExt2Int[familyJ], familiesExt2Int[familyI]] = familySimilarityIJ

                # if we have the identity matrix here, it means all GCFs are separate
                # (nothing to cluster). Note: still can crash if values are
                # sufficiently low. Catch this error later
                #if np.count_nonzero(simMatrix) == 0:
                    #clanLabels = []
                    #continue

                # add main diagonal
                for family in range(len(familyIdx)):
                    famSimMatrix[family,family] = 1.0

                bs_similarity_families = famSimMatrix.tolist()

                #clanLabels = pysapc.SAP(damping=damping, max_iter=500,
                                    #preference='min').fit_predict(famSimMatrix)
                af = AffinityPropagation(damping=damping, max_iter=1000, convergence_iter=200, affinity="precomputed").fit(famSimMatrix)
                labelsClans = af.labels_
                exemplarsClans = af.cluster_centers_indices_

                # affinity propagation can fail in some circumstances (e.g. only singletons)
                if exemplarsClans is not None:
                    # translate and record GCF label instead of GCF number
                    clanLabels = [familyIdx[exemplarsClans[labelsClans[i]]] for i in range(len(familyIdx))]
                else:
                    clanLabels = []

        else:
            clanLabels = []

        if len(clanLabels) > 0:
            clansDict = defaultdict(list)
            for i in range(len(familyIdx)):
                clansDict[clanLabels[i]].append(familyIdx[i])

            fam2clan = dict(zip(familyIdx,clanLabels))

            bs_families = [{"id": "FAM_{:05d}".format(family),
                            'members': [bgcExt2Int[member] for member in members],
                            "clan": "CLAN_{:03d}".format(fam2clan[family])}
                           for family, members in familiesDict.items()]
            bs_clans = [{"id": "CLAN_{:03d}".format(clan), 'members': members}
                           for clan, members in clansDict.items()]
        else:
            bs_families = [{"id": "FAM_{:05d}".format(family),
                            'members': [bgcExt2Int[member] for member in members], }
                           for family, members in familiesDict.items()]

        family_data["families"] = []
        for family, members in familiesDict.items():
            family_data["families"].append({
                "label": "FAM_{:05d}".format(family),
                "members": members # use external indexing
            })

        # Positional alignment information is based on DomainCountGene, which
        # does not contain empty genes (i.e. with no domains).
        domainGenes2allGenes = {}

        ## BGC Family alignment information
        bs_families_alignment = []
        for family, members in familiesDict.items():
            for bgc in members:
                domainGenes2allGenes[bgc] = {}
                has_domains = 0
                for orf in range(len(bs_data[bgcExt2Int[bgc]]["orfs"])):
                    if len(bs_data[bgcExt2Int[bgc]]["orfs"][orf]["domains"]) > 0:
                        domainGenes2allGenes[bgc][has_domains] = orf
                        has_domains += 1

            assert (len(members) > 0), "Error: bs_families[{}] have no members, something went wrong?".format(family)

            ref_genes_ = set()
            aln = []


            for bgc in members:
                if bgc == family:
                    aln.append([ [gene_num, 0] for gene_num in range(len(bs_data[bgcExt2Int[family]]["orfs"]))])
                else:
                    try:
                        a, b, length, reverse = pos_alignments[family][bgc]
                    except:
                        b, a, length, reverse = pos_alignments[bgc][family]

                        if length == 0:
                            pass
                        elif reverse:
                            # special case. bgc was reference (first) in lcs
                            a = domainGenes2allGenes[family][len(bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[family]].gene_domain_counts)-a-length]
                            b = domainGenes2allGenes[bgc][b+length-1] # -1 go to 0-index
                        else:
                            a = domainGenes2allGenes[family][a]
                            b = domainGenes2allGenes[bgc][b]
                    else:
                        a = domainGenes2allGenes[family][a]
                        if length == 0:
                            pass

                        elif reverse:

                            b = domainGenes2allGenes[bgc][len(bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].gene_domain_counts)-b-1]
                        else:
                            b = domainGenes2allGenes[bgc][b]


                    if length == 0:
                        length = 1
                        # let's try aligning using the genes with most domains
                        # after all, they ended up being in the same GCF
                        # for some reason
                        x = max(bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[family]].gene_domain_counts)
                        x = bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[family]].gene_domain_counts.index(x)
                        a = domainGenes2allGenes[family][x]

                        y = max(list(bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].gene_domain_counts))
                        y = bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].gene_domain_counts.index(y)

                        #check orientation
                        if bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[family]].gene_orientations[x] == bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].gene_orientations[y]:
                            b = domainGenes2allGenes[bgc][y]
                            reverse = False
                        else:
                            b = domainGenes2allGenes[bgc][len(bgc_collection.bgc_collection_dict[bgc_collection.bgc_name_tuple[bgc]].gene_domain_counts)-y-1]
                            reverse = True

                    ref_genes_.add(a)
                    bgc_algn = []

                    for gene_num in range(len(bs_data[bgcExt2Int[bgc]]["orfs"])):
                        if gene_num == b: # this is the reference gene for this bgc
                            if reverse:
                                bgc_algn.append([a, -100])
                            else:
                                bgc_algn.append([a, 100])
                        else:
                            bgc_algn.append([-1, 100])

                    aln.append(bgc_algn)

            ref_genes = list(ref_genes_)

            fam_alignment = {
                "id" : "FAM_{:05d}".format(family),
                "ref" : bgcExt2Int[family],
                "newick" : newick_trees[family],
                "ref_genes" : ref_genes,
                "aln" : aln
            }
            bs_families_alignment.append(fam_alignment)
        ## End of BGC Family alignment information

        # column1: BGC, column2: clustering pseudo family
        
        logging.debug("  Writing clustering file")
        clustering_file_path = os.path.join(pathBase, "{}_clustering_c{:4.2f}.tsv".format(className, cutoff))
        with open(clustering_file_path, "w") as clustering_file:
            clustering_file.write('#BGC Name\tFamily Number\n')
            for family in familyIdx:
                for bgc in familiesDict[family]:
                    clustering_file.write("{}\t{}\n".format(bgc_collection.bgc_name_tuple[bgc], family))

        # get family-family similarity matrix
        bs_similarity_families = [[get_composite_bgc_similarities([bgcs[bid] for bid in bs_families[row]["members"]], [bgcs[bid] for bid in bs_families[col]["members"]], simDict) if (row != col) else (1.00, (1.00, bgcs[bs_families[row]["members"][0]], bgcs[bs_families[row]["members"][0]]), (1.00, bgcs[bs_families[row]["members"][0]], bgcs[bs_families[row]["members"][0]])) for col in
                         range(row+1)] for row in range(len(bs_families))]

        family_data["families_similarity"] = bs_similarity_families

        ## Write html output folder structure (and update bigscape_results.js) for this module
        htmlFolder_run = "{}_c{:.2f}".format(htmlFolder, cutoff)
        assert os.path.isdir(htmlFolder_run)
        module_html_path = os.path.join(htmlFolder_run, className)
        create_directory(module_html_path, "Network HTML", False)
        with open(os.path.join(module_html_path, "bs_data.js"), "w") as bs_data_js:
            bs_data_js.write("var bs_data={};\n".format(json.dumps(bs_data, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_data_js.write("dataLoaded('bs_data');\n")
            # TODO: better file path handling
        shutil.copy(os.path.join(os.path.realpath("./"), "html_template", "index_html"), os.path.join(module_html_path, "index.html"))

        ## Write bgc_networks.js
        with open(os.path.join(module_html_path, "bs_networks.js"), "w") as bs_networks_js:
            bs_networks_js.write("var bs_similarity={};\n".format(json.dumps(bs_distances, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families={};\n".format(json.dumps(bs_families, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_families_alignment={};\n".format(json.dumps(bs_families_alignment, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("var bs_similarity_families={};\n".format(json.dumps(bs_similarity_families, indent=4, separators=(',', ':'), sort_keys=True)))
            if len(clanLabels) > 0:
                bs_networks_js.write("var bs_clans={};\n".format(json.dumps(bs_clans, indent=4, separators=(',', ':'), sort_keys=True)))
            bs_networks_js.write("dataLoaded('bs_networks');\n")


        if len(clanLabels) > 0:
            logging.debug("   Writing Clans file")
            clans_file_path = os.path.join(pathBase, "{}_clans_{:4.2f}_{:4.2f}.tsv".format(className,clanClassificationCutoff,clanDistanceCutoff))
            with open(clans_file_path,'w') as clansFile:
                clansFile.write('#BGC Name\tClan Number\tFamily Number\n')
                for clan in clansDict.keys():
                    for family in clansDict[clan]:
                        for bgc in familiesDict[family]:
                            clansFile.write("{}\t{}\t{}\n".format(bgc_collection.bgc_name_tuple[bgc], clan, family))

        results[htmlFolder_run] = family_data

    return results
