import os
import numpy as np

from collections import defaultdict
from difflib import SequenceMatcher
from Bio import pairwise2
from Bio.SubsMat.MatrixInfo import pam250 as scoring_matrix
from scipy.optimize import linear_sum_assignment


from src.big_scape.bgc_info import BgcInfo
from src.utility.fasta import fasta_parser

def gen_unrelated_pair_distance(run, cluster_a: BgcInfo, cluster_b: BgcInfo):

    num_non_anchor_domains = 0
    num_anchor_domains = 0
    # Count total number of anchor and non-anchor domain to report in the
    # network file. Apart from that, these BGCs are totally unrelated.
    for domain in cluster_a.ordered_domain_set:
        # This is a bit of a hack. If pfam domain ids ever change in size
        # we'd be in trouble. The previous approach was to .split(".")[0]
        # but it's more costly
        if domain[:7] in run.network.anchor_domains:
            num_anchor_domains += len(cluster_a.domain_name_info[domain])
        else:
            num_non_anchor_domains += len(cluster_a.domain_name_info[domain])

    for domain in cluster_b.ordered_domain_set:
        if domain[:7] in run.network.anchor_domains:
            num_anchor_domains += len(cluster_b.domain_name_info[domain])
        else:
            num_non_anchor_domains += len(cluster_b.domain_name_info[domain])

    return 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, num_non_anchor_domains, num_anchor_domains, 0, 0, 0, 0

def get_lcs(a_string, b_string, a_match_len, b_match_len):
    # Longest Common Substring (LCS)
    # construct object for finding LCS
    seqmatch = SequenceMatcher(None, a_string, b_string)
    # find the LCS
    match = seqmatch.find_longest_match(0, a_match_len, 0, b_match_len)
    # unpack
    a_start, b_start, match_length = match
    return a_start, b_start, match_length

def get_lcs_fwd(cluster_a, cluster_b):
    return get_lcs(cluster_a.gene_string, cluster_b.gene_string, cluster_a.num_genes, cluster_b.num_genes)

def get_lcs_rev(cluster_a, cluster_b):
    return get_lcs(cluster_a.gene_string, cluster_b.gene_string_rev, cluster_a.num_genes, cluster_b.num_genes)

def score_expansion(x_string_, y_string_, downstream):
    """
    Input:
    x_string: list of strings. This one tries to expand based on max score
    y_string: reference list. This was already expanded.
    downstream: If true, expansion goes from left to right.

    Output:
    max_score, a
    Where a is length of expansion.
    """
    match = 5
    mismatch = -3
    gap = -2

    if downstream:
        x_string = x_string_
        y_string = y_string_
    else:
        # Expansion goes upstream. It's easier to flip both slices and proceed
        # as if going downstream.
        x_string = list(reversed(x_string_))
        y_string = list(reversed(y_string_))


    # how many gaps to open before calling a mismatch is more convenient?
    # max_gaps_before_mismatch = abs(int(match/gap))
    max_score = 0
    score = 0

    pos_y = 0
    expand_len = 0

    for pos_x, gene in enumerate(x_string):
        # try to find g within the rest of the slice
        # This has the obvious problem of what to do if a gene is found _before_
        # the current y_slice (translocation). Duplications could also complicate
        # things
        try:
            match_pos = y_string[pos_y:].index(gene)
        except ValueError:
            score += mismatch
        else:
            score += match + match_pos*gap
            pos_y += match_pos + 1 # move pointer one position after match

            # Greater or equals. 'equals' because even if max_score isn't
            # larger, it's an expansion
            if score >= max_score:
                max_score = score
                # keep track of expansion. Account for zero-based numbering
                expand_len = pos_x + 1

        #print(pos_x, score, status, g)
        #print("")

    return max_score, expand_len

def calc_jaccard(intersect, overlap):
    return len(intersect) / len(overlap)

def calc_adj_idx(cluster_a, cluster_b,
                 cluster_a_dom_start, cluster_b_dom_start,
                 cluster_a_dom_end, cluster_b_dom_end):
    # ADJACENCY INDEX
    # calculates the Tanimoto similarity of pairs of adjacent domains

    if len(cluster_a.ordered_domain_list[cluster_a_dom_start:cluster_a_dom_end]) < 2 or len(cluster_b.ordered_domain_list[cluster_b_dom_start:cluster_b_dom_end]) < 2:
        adj_idx = 0.0
    else:
        set_a_pairs = set()
        for cluster_idx in range(cluster_a_dom_start, cluster_a_dom_end-1):
            set_a_pairs.add(tuple(sorted([cluster_a.ordered_domain_list[cluster_idx], cluster_a.ordered_domain_list[cluster_idx+1]])))

        set_b_pairs = set()
        for cluster_idx in range(cluster_b_dom_start, cluster_b_dom_end-1):
            set_b_pairs.add(tuple(sorted([cluster_b.ordered_domain_list[cluster_idx], cluster_b.ordered_domain_list[cluster_idx+1]])))

        # same treatment as in Jaccard
        intersect = set_a_pairs & set_b_pairs
        overlap = set_a_pairs | set_b_pairs
        adj_idx = calc_jaccard(intersect, overlap)
    return adj_idx

def process_orientation(cluster_a, cluster_b):
    # get LCS for this pair
    a_start, b_start, match_length = get_lcs_fwd(cluster_a, cluster_b)

    a_start_rev, b_start_rev, match_length_rev = get_lcs_rev(cluster_a, cluster_b)

    # We need to keep working with the correct orientation
    # forward case
    if match_length > match_length_rev or (match_length == match_length_rev and cluster_a.gene_orientations[a_start] == cluster_b.gene_orientations[b_start]):
        use_b_string = cluster_b.gene_string
        # note: these slices are in terms of genes, not domains (which are
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        slice_start_a = a_start
        slice_start_b = b_start
        slice_length_a = match_length
        slice_length_b = match_length

        reverse = False
    # reverse based on match length case
    elif match_length < match_length_rev:
        use_b_string = cluster_b.gene_string_rev

        slice_start_a = a_start_rev
        slice_start_b = b_start_rev
        slice_length_a = match_length_rev
        slice_length_b = match_length_rev

        # We'll need to know if we're working in reverse so that the start
        # postion of the final slice can be transformed to the original orientation
        reverse = True

    # only one gene matches case
    # choose the one with the most domains
    elif match_length == 1:
        seqmatch = SequenceMatcher(None, cluster_a.gene_string, cluster_b.gene_string)
        max_domains = 0
        index_a = 0   # index in A with the gene with most domains
        index_b = 0   # index in B with the gene with most domains
        for a_start, b_start, block_match_len in seqmatch.get_matching_blocks():
            if block_match_len != 0:
                if cluster_a.gene_domain_counts[a_start] > max_domains:
                    index_a = a_start
                    index_b = b_start
                    max_domains = cluster_a.gene_domain_counts[a_start]


        # note: these slices are in terms of genes, not domains (which are
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        slice_start_a = index_a
        slice_length_a = 1
        slice_length_b = 1

        if cluster_a.gene_orientations[index_a] == cluster_b.gene_orientations[index_b]:
            slice_start_b = index_b
            use_b_string = cluster_b.gene_string
            reverse = False
        else:
            slice_start_b = len(cluster_b.gene_orientations) - index_b - 1
            use_b_string = cluster_b.gene_string_rev
            reverse = True

    # s == sr and (s == 0 or s > 1)
    # match_len == match_len_rev, and:
    # match_len == 0, or match_len > 1
    else:
        use_b_string = cluster_b.gene_string
        # note: these slices are in terms of genes, not domains (which are
        # ultimately what is used for distance)
        # Currently, the following values represent the Core Overlap
        slice_start_a = a_start
        slice_start_b = b_start
        slice_length_a = match_length
        slice_length_b = match_length

        reverse = False
    return slice_start_a, slice_start_b, slice_length_a, slice_length_b, use_b_string, reverse

def calc_distance_lcs(run, cluster_a: BgcInfo, cluster_b: BgcInfo, weights,
                      aligned_domain_sequences):
    """Compare two clusters using information on their domains, and the
    sequences of the domains.
    This version first tries to search for the largest common slices of both BGCs by
    finding the Longest Common Substring of genes (based on domain content), also
    trying the reverse on one of the BGCs (capital "B" will be used when the 'true'
    orientation of the BGC is found)
    Then the algorithm will try to expand the slices to either side. Once each
    slice are found, they have to pass one more validation: a core biosynthetic
    gene (marked in antiSMASH) must be present.
    Capital letters indicate the "true" orientation of the BGC (first BGC, 'A'
    is kept fixed)

    Output:
    raw distance, jaccard index, domain sequence similarity, adjacency index,
    raw DSS non-anchor, raw DSS anchor, number of non-anchor domains in the pair
    (S), number of anchor domains in the pair

    """

    # FIXME: some variables were renamed in such a way that the funcionality of this function
    # is probably broken. See line 163 for example. this used to be an assignment to a different
    # function, now it's the same one. Probably because uppercase characters were replaced with
    # lowercase.
    # uppercase characters may be the reverse version of whatever BGC is being referred to

    ## temporary cluster variables
    cluster_a_dom_start = 0
    cluster_a_dom_end = len(cluster_a.ordered_domain_list)

    cluster_b_dom_start = 0
    cluster_b_dom_end = len(cluster_b.ordered_domain_list)

    # initialize domain sequence slices
    # They might change if we manage to find a valid overlap
    cluster_a_domain_seq_slice_bottom = defaultdict(int)
    cluster_a_domain_seq_slice_top = defaultdict(int)
    for domain in cluster_a.ordered_domain_set:
        cluster_a_domain_seq_slice_bottom[domain] = 0
        cluster_a_domain_seq_slice_top[domain] = len(cluster_a.domain_name_info[domain])

    cluster_b_domain_seq_slice_bottom = defaultdict(int)
    cluster_b_domain_seq_slice_top = defaultdict(int)
    for domain in cluster_b.ordered_domain_set:
        cluster_b_domain_seq_slice_bottom[domain] = 0
        cluster_b_domain_seq_slice_top[domain] = len(cluster_b.domain_name_info[domain])

    cluster_b_reverse_gene_domain_counts = list(reversed(cluster_b.gene_domain_counts))

    cluster_a_temp_domain_set = cluster_a.ordered_domain_set
    cluster_b_temp_domain_set = cluster_b.ordered_domain_set

    ## common variables
    # unpack weights
    jaccard_weight, dss_weight, ai_weight, anchor_boost = weights

    temp_domain_fastas = {}

    intersect = cluster_a.ordered_domain_set & cluster_b.ordered_domain_set

    # Detect totally unrelated pairs from the beginning
    # this returns early if there are unrelated pairs
    if len(intersect) == 0:
        return gen_unrelated_pair_distance(run, cluster_a, cluster_b)


    slice_data = process_orientation(cluster_a, cluster_b)

    slice_start_a, slice_start_b, slice_length_a, slice_length_b, use_b_string, reverse = slice_data


    lcs_start_a = slice_start_a
    lcs_start_b = slice_start_b
    seed_length = slice_length_a

    is_glocal = run.options.mode == "glocal"
    is_auto = run.options.mode == "auto"
    is_contig_edge_a = cluster_a.bgc_info.contig_edge
    is_contig_edge_b = cluster_b.bgc_info.contig_edge

    if is_glocal or (is_auto and (is_contig_edge_a or is_contig_edge_b)):
        #X: bgc that drive expansion
        #Y: the other bgc
        # forward: True if expansion is to the right
        # returns max_score, final positions for X and Y
        #score_expansion(X_string, dcg_X, Y_string, dcg_Y, downstream=True/False)

        # Expansion is relatively costly. We ask for a minimum of 3 genes
        # for the core overlap before proceeding with expansion.
        biosynthetic_hit_a = False
        for biosynthetic_position in cluster_a.bio_synth_core_positions:
            if biosynthetic_position >= slice_start_a and biosynthetic_position <= (slice_start_a+slice_length_a):
                biosynthetic_hit_a = True
                break

        if slice_length_a >= 3 or biosynthetic_hit_a:
            # LEFT SIDE
            # Find which bgc has the least genes to the left. If both have the same
            # number, find the one that drives the expansion with highest possible score
            if slice_start_a == slice_start_b:
                # assume complete expansion of A, try to expand B
                exp_score_b, expansion_len_b = score_expansion(use_b_string[:slice_start_b], cluster_a.gene_string[:slice_start_a], False)
                # assume complete expansion of B, try to expand A
                exp_score_a, expansion_len_a = score_expansion(cluster_a.gene_string[:slice_start_a], use_b_string[:slice_start_b], False)

                if exp_score_a > exp_score_b or (exp_score_a == exp_score_b and expansion_len_a > expansion_len_b):
                    slice_length_a += expansion_len_a
                    slice_length_b += len(use_b_string[:slice_start_b])
                    slice_start_a -= expansion_len_a
                    slice_start_b = 0
                else:
                    slice_length_a += len(cluster_a.gene_string[:slice_start_a])
                    slice_length_b += expansion_len_b
                    slice_start_a = 0
                    slice_start_b -= expansion_len_b

            else:
                # A is shorter upstream. Assume complete extension. Find B's extension
                if slice_start_a < slice_start_b:
                    exp_score_b, exp_length_b = score_expansion(use_b_string[:slice_start_b], cluster_a.gene_string[:slice_start_a], False)

                    slice_length_a += len(cluster_a.gene_string[:slice_start_a])
                    slice_length_b += exp_length_b
                    slice_start_a = 0
                    slice_start_b -= exp_length_b
                else:
                    exp_score_a, exp_length_a = score_expansion(cluster_a.gene_string[:slice_start_a], use_b_string[:slice_start_b], False)

                    slice_length_a += exp_length_a
                    slice_length_b += len(use_b_string[:slice_start_b])
                    slice_start_a -= exp_length_a
                    slice_start_b = 0

            # RIGHT SIDE
            # check which side is the shortest downstream. If both BGCs have the same
            # length left, choose the one with the best expansion
            downstream_a = cluster_a.num_genes - slice_start_a - slice_length_a
            downstream_b = cluster_b.num_genes - slice_start_b - slice_length_b
            if downstream_a == downstream_b:
                # assume complete extension of A, try to expand B
                exp_score_b, exp_len_b = score_expansion(use_b_string[slice_start_b+slice_length_b:], cluster_a.gene_string[slice_start_a+slice_length_a:], True)
                # assume complete extension of B, try to expand A
                exp_score_a, exp_len_a = score_expansion(cluster_a.gene_string[slice_start_a+slice_length_a:], use_b_string[slice_start_b+slice_length_b:], True)


                if (exp_score_a == exp_score_b and exp_len_a > exp_len_b) or exp_score_a > exp_score_b:
                    slice_length_a += exp_len_a
                    slice_length_b += len(use_b_string[slice_start_b+slice_length_b:])
                else:
                    slice_length_a += len(cluster_a.gene_string[slice_start_a+slice_length_a:])
                    slice_length_b += exp_len_b

            else:
                if downstream_a < downstream_b:
                    # extend all of remaining A
                    exp_score_b, exp_len_b = score_expansion(use_b_string[slice_start_b+slice_length_b:], cluster_a.gene_string[slice_start_a+slice_length_a:], True)

                    slice_length_a += len(cluster_a.gene_string[slice_start_a+slice_length_a:])
                    slice_length_b += exp_len_b

                else:
                    exp_score_a, exp_len_a = score_expansion(cluster_a.gene_string[slice_start_a+slice_length_a:], use_b_string[slice_start_b+slice_length_b:], True)
                    slice_length_a += exp_len_a
                    slice_length_b += len(use_b_string[slice_start_b+slice_length_b:])

        if min(slice_length_a, slice_length_b) >= 5:
            # First test passed. Find if there is a biosynthetic gene in both slices
            # (note that even if they are, currently we don't check whether it's
            # actually the _same_ gene)
            biosynthetic_hit_a = False
            biosynthetic_hit_b = False

            for biosynthetic_position in cluster_a.bio_synth_core_positions:
                if biosynthetic_position >= slice_start_a and biosynthetic_position <= (slice_start_a+slice_length_a):
                    biosynthetic_hit_a = True
                    break


            # return to original orientation if needed
            if reverse:
                slice_start_b = cluster_b.num_genes - slice_start_b - slice_length_b

            # using original bio_synth_core_positions
            for biosynthetic_position in cluster_b.bio_synth_core_positions:
                if biosynthetic_position >= slice_start_b and biosynthetic_position <= (slice_start_b + slice_length_b):
                    biosynthetic_hit_b = True
                    break

            # finally...
            if biosynthetic_hit_a and biosynthetic_hit_b:
                cluster_a_dom_start = sum(cluster_a.gene_domain_counts[:slice_start_a])
                cluster_a_dom_end = cluster_a_dom_start + sum(cluster_a.gene_domain_counts[slice_start_a:slice_start_a+slice_length_a])
                cluster_a_temp_domain_set = set(cluster_a.ordered_domain_list[cluster_a_dom_start:cluster_a_dom_end])

                if reverse:
                    cluster_b_dom_start = sum(cluster_b_reverse_gene_domain_counts[:slice_start_b])
                    cluster_b_dom_end = cluster_b_dom_start + sum(cluster_b_reverse_gene_domain_counts[slice_start_b:slice_start_b+slice_length_b])
                else:
                    cluster_b_dom_start = sum(cluster_b.domain_gene_counts[:slice_start_b])
                    cluster_b_dom_end = cluster_b_dom_start + sum(cluster_b.domain_gene_counts[slice_start_b:slice_start_b+slice_length_b])
                cluster_b_temp_domain_set = set(cluster_b.ordered_domain_list[cluster_b_dom_start:cluster_b_dom_end])

                intersect = cluster_a_temp_domain_set & cluster_b_temp_domain_set

                # re-adjust the indices for each domain so we get only the sequence
                # tags in the selected slice. First step: find out which is the
                # first copy of each domain we're using
                for domain in cluster_a.ordered_domain_list[:cluster_a_dom_start]:
                    cluster_a_domain_seq_slice_bottom[domain] += 1
                for domain in cluster_b.ordered_domain_list[:cluster_b_dom_start]:
                    cluster_b_domain_seq_slice_bottom[domain] += 1

                # Step 2: work with the last copy of each domain.
                # Step 2a: make top = bottom
                for domain in cluster_a_temp_domain_set:
                    cluster_a_domain_seq_slice_top[domain] = cluster_a_domain_seq_slice_bottom[domain]
                for domain in cluster_b_temp_domain_set:
                    cluster_b_domain_seq_slice_top[domain] = cluster_b_domain_seq_slice_bottom[domain]

                # Step 2b: increase top with the domains in the slice
                for domain in cluster_a.ordered_domain_list[cluster_a_dom_start:cluster_a_dom_end]:
                    cluster_a_domain_seq_slice_top[domain] += 1
                for domain in cluster_b.ordered_domain_list[cluster_b_dom_start:cluster_b_dom_end]:
                    cluster_b_domain_seq_slice_top[domain] += 1

            #else:
                #print(" - - Not a valid overlap found - - (no biosynthetic genes)\n")

        #else:
                #print(" - - Not a valid overlap found - - (shortest slice not large enough)\n")


    # JACCARD INDEX
    union = cluster_a.ordered_domain_set | cluster_b.ordered_domain_set
    jaccard_index = calc_jaccard(intersect, union)


    # DSS INDEX
    #domain_difference: Difference in sequence per domain. If one cluster does
    # not have a domain at all, but the other does, this is a (complete)
    # difference in sequence 1. If both clusters contain the domain once, and
    # the sequence is the same, there is a seq diff of 0.
    #S: Max occurence of each domain
    domain_difference_anchor, num_anchor_domains = 0, 0
    domain_difference, num_non_anchor_domains = 0, 0


    not_intersect = cluster_a_temp_domain_set.symmetric_difference(cluster_b_temp_domain_set)

    # Case 1
    #no need to look at seq identity, since these domains are unshared
    for unshared_domain in not_intersect:
        #for each occurence of an unshared domain do domain_difference += count
        # of domain and S += count of domain

        try:
            num_unshared = cluster_a_domain_seq_slice_top[unshared_domain] - cluster_a_domain_seq_slice_bottom[unshared_domain]
        except KeyError:
            num_unshared = cluster_b_domain_seq_slice_top[unshared_domain] - cluster_b_domain_seq_slice_bottom[unshared_domain]

        # don't look at domain version, hence the split
        if unshared_domain[:7] in run.network.anchor_domains:
            domain_difference_anchor += num_unshared
        else:
            domain_difference += num_unshared

    num_non_anchor_domains = domain_difference # can be done because it's the first use of these
    num_anchor_domains = domain_difference_anchor

    # Cases 2 and 3 (now merged)
    missing_aligned_domain_files = []

    for shared_domain in intersect:
        specific_domain_list_a = cluster_a.domain_name_info[shared_domain]
        specific_domain_list_b = cluster_b.domain_name_info[shared_domain]

        num_copies_a = cluster_a_domain_seq_slice_top[shared_domain] - cluster_a_domain_seq_slice_bottom[shared_domain]
        num_copies_b = cluster_b_domain_seq_slice_top[shared_domain] - cluster_b_domain_seq_slice_bottom[shared_domain]

        temp_domain_fastas.clear()

        accumulated_distance = 0

        # Fill distance matrix between domain's A and B versions
        distance_matrix = np.ndarray((num_copies_a, num_copies_b))
        for domsa in range(num_copies_a):
            for domsb in range(num_copies_b):
                sequence_tag_a = specific_domain_list_a[domsa + cluster_a_domain_seq_slice_bottom[shared_domain]]
                sequence_tag_b = specific_domain_list_b[domsb + cluster_b_domain_seq_slice_bottom[shared_domain]]

                seq_length = 0
                matches = 0
                gaps = 0

                try:
                    aligned_seq_a = aligned_domain_sequences[sequence_tag_a]
                    aligned_seq_b = aligned_domain_sequences[sequence_tag_b]

                except KeyError:
                    # For some reason we don't have the multiple alignment files.
                    # Try manual alignment
                    if shared_domain not in missing_aligned_domain_files and run.options.verbose:
                        # this will print everytime an unfound <domain>.algn is not found for every
                        # distance calculation (but at least, not for every domain pair!)
                        print("  Warning: {}.algn not found. Trying pairwise alignment...".format(shared_domain))
                        missing_aligned_domain_files.append(shared_domain)

                    try:
                        unaligned_seq_a = temp_domain_fastas[sequence_tag_a]
                        unaligned_seq_b = temp_domain_fastas[sequence_tag_b]
                    except KeyError:
                        # parse the file for the first time and load all the sequences
                        with open(os.path.join(run.directories.domains, shared_domain + ".fasta"),"r") as domain_fasta_handle:
                            temp_domain_fastas = fasta_parser(domain_fasta_handle)

                        unaligned_seq_a = temp_domain_fastas[sequence_tag_a]
                        unaligned_seq_b = temp_domain_fastas[sequence_tag_b]

                    # gap_open = -15
                    # gap_extend = -6.67. These parameters were set up by Emzo
                    align_score = pairwise2.align.globalds(unaligned_seq_a, unaligned_seq_b, scoring_matrix, -15, -6.67, one_alignment_only=True)
                    best_alignment = align_score[0]
                    aligned_seq_a = best_alignment[0]
                    aligned_seq_b = best_alignment[1]


                # - Calculate aligned domain sequences similarity -
                # Sequences *should* be of the same length unless something went
                # wrong elsewhere
                if len(aligned_seq_a) != len(aligned_seq_b):
                    print("\tWARNING: mismatch in sequences' lengths while calculating sequence identity ({})".format(shared_domain))
                    print("\t  Specific domain 1: {} len: {}".format(sequence_tag_a, str(len(aligned_seq_a))))
                    print("\t  Specific domain 2: {} len: {}".format(sequence_tag_b, str(len(aligned_seq_b))))
                    seq_length = min(len(aligned_seq_a), len(aligned_seq_b))
                else:
                    seq_length = len(aligned_seq_a)

                for position in range(seq_length):
                    if aligned_seq_a[position] == aligned_seq_b[position]:
                        if aligned_seq_a[position] != "-":
                            matches += 1
                        else:
                            gaps += 1

                distance_matrix[domsa][domsb] = 1 - (matches/(seq_length-gaps))

        #Only use the best scoring pairs
        best_indexes = linear_sum_assignment(distance_matrix)
        accumulated_distance = distance_matrix[best_indexes].sum()

        # the difference in number of domains accounts for the "lost" (or not duplicated) domains
        sum_seq_dist = (abs(num_copies_a-num_copies_b) + accumulated_distance)  #essentially 1-sim
        normalization_element = max(num_copies_a, num_copies_b)

        if shared_domain[:7] in run.network.anchor_domains:
            num_anchor_domains += normalization_element
            domain_difference_anchor += sum_seq_dist
        else:
            num_non_anchor_domains += normalization_element
            domain_difference += sum_seq_dist

    if num_anchor_domains != 0 and num_non_anchor_domains != 0:
        dss_non_anchor = domain_difference / num_non_anchor_domains
        dss_anchor = domain_difference_anchor / num_anchor_domains

        # Calculate proper, proportional weight to each kind of domain
        non_anchor_prct = num_non_anchor_domains / (num_non_anchor_domains + num_anchor_domains)
        anchor_prct = num_anchor_domains / (num_non_anchor_domains + num_anchor_domains)

        # boost anchor subcomponent and re-normalize
        non_anchor_weight = non_anchor_prct / (anchor_prct*anchor_boost + non_anchor_prct)
        anchor_weight = anchor_prct*anchor_boost / (anchor_prct*anchor_boost + non_anchor_prct)

        # Use anchorboost parameter to boost percieved rDSS_anchor
        dss = (non_anchor_weight*dss_non_anchor) + (anchor_weight*dss_anchor)

    elif num_anchor_domains == 0:
        dss_non_anchor = domain_difference / num_non_anchor_domains
        dss_anchor = 0.0

        dss = dss_non_anchor

    else: #only anchor domains were found
        dss_non_anchor = 0.0
        dss_anchor = domain_difference_anchor / num_anchor_domains

        dss = dss_anchor

    dss = 1-dss #transform into similarity

    adj_index = calc_adj_idx(cluster_a, cluster_b,
                             cluster_a_dom_start, cluster_b_dom_start,
                             cluster_a_dom_end, cluster_b_dom_end)


    distance = 1 - (jaccard_weight * jaccard_index) - (dss_weight * dss) - (ai_weight * adj_index)

    # This could happen due to numerical innacuracies
    if distance < 0.0:
        if distance < -0.000001: # this definitely is something else...
            print("Negative distance detected!")
            print(distance)
            print("{} - {}".format(cluster_a.name, cluster_b.name))
            print("J: {}\tDSS: {}\tAI: {}".format(str(jaccard_index), str(dss), str(adj_index)))
            print("Jw: {}\tDSSw: {}\tAIw: {}".format(str(jaccard_weight), str(dss_weight), str(ai_weight)))
        distance = 0.0

    rev = 0.0
    if reverse:
        rev = 1.0

    return distance, jaccard_index, dss, adj_index, dss_non_anchor, dss_anchor, num_non_anchor_domains, num_anchor_domains, lcs_start_a, lcs_start_b, seed_length, rev
