"""Module to contain the BgcDomainInfo class, used in distance calculation
between two pairs specifically

Authors: Jorge Navarro, Arjan Draisma
"""

from collections import defaultdict

from src.big_scape.scores import score_expansion
from src.big_scape.bgc_info import BgcInfo


class BgcDomainInfo():
    """Class to contain bgc specific domain information used in distance
    calculation
    """
    def __init__(self, bgc_a: BgcInfo, bgc_b: BgcInfo):

        self.intersect = bgc_a.ordered_domain_set & bgc_b.ordered_domain_set

        self.a_dom_start = 0
        self.a_dom_end = len(bgc_a.ordered_domain_list)

        self.b_dom_start = 0
        self.b_dom_end = len(bgc_b.ordered_domain_list)

        # initialize domain sequence slices
        # They might change if we manage to find a valid overlap
        self.a_dom_seq_slice_bot = defaultdict(int)
        self.a_dom_seq_slice_top = defaultdict(int)
        for domain in bgc_a.ordered_domain_set:
            self.a_dom_seq_slice_bot[domain] = 0
            self.a_dom_seq_slice_top[domain] = len(bgc_a.domain_name_info[domain])

        self.b_dom_seq_slice_bot = defaultdict(int)
        self.b_dom_seq_slice_top = defaultdict(int)
        for domain in bgc_b.ordered_domain_set:
            self.b_dom_seq_slice_bot[domain] = 0
            self.b_dom_seq_slice_top[domain] = len(bgc_b.domain_name_info[domain])

        self.b_rev_gene_dom_counts = list(reversed(bgc_b.gene_domain_counts))

        self.a_dom_set = bgc_a.ordered_domain_set
        self.b_dom_set = bgc_b.ordered_domain_set


    def expand_score(self, run, bgc_a, bgc_b, slice_data):
        """Expand the selection of genes for score calculation. Since this
        function is relatively costly, Score is only expanded if LCS is at
        least 3 genes long, or has at least one gene marked as core
        biosynthetic

        Inputs:
            run: run details for this execution of BiG-SCAPE
            cluster_a, cluster_b: bgcInfo objects for bgcs of which distance is
            calculated
            slice_data: slice data from the bgc comparison
        """
        slice_start_a = slice_data[0]
        slice_start_b = slice_data[1]
        slice_length_a = slice_data[2]
        slice_length_b = slice_data[3]
        use_b_string = slice_data[4]
        reverse = slice_data[5]

        is_glocal = run.options.mode == "glocal"
        is_auto = run.options.mode == "auto"

        is_contig_edge_a = bgc_a.bgc_data.contig_edge
        is_contig_edge_b = bgc_b.bgc_data.contig_edge
        is_contig_edge = is_contig_edge_a or is_contig_edge_b

        # don't do anything if we're not in glocal mode or we're not a contig
        # edge
        if not is_glocal and not (is_auto and is_contig_edge):
            return

        # Expansion is relatively costly. We ask for a minimum of 3 genes
        # for the core overlap before proceeding with expansion.
        biosynthetic_hit_a = False
        for biosynthetic_position in bgc_a.bio_synth_core_positions:
            after_start = biosynthetic_position >= slice_start_a

            slice_end_a = slice_start_a+slice_length_a
            before_end = biosynthetic_position <= slice_end_a

            # only count biosynthetic hits within the slice range
            if after_start and before_end:
                biosynthetic_hit_a = True
                break

        if slice_length_a >= 3 or biosynthetic_hit_a:
            # LEFT SIDE
            # Find which bgc has the least genes to the left. If both have the
            # same number, find the one that drives the expansion with highest
            # possible score
            if slice_start_a == slice_start_b:
                # assume complete expansion of A, try to expand B
                exp_score_b, expansion_len_b = score_expansion(
                    use_b_string[:slice_start_b],
                    bgc_a.gene_string[:slice_start_a],
                    False
                )
                # assume complete expansion of B, try to expand A
                exp_score_a, expansion_len_a = score_expansion(
                    bgc_a.gene_string[:slice_start_a],
                    use_b_string[:slice_start_b],
                    False
                )

                if (
                    exp_score_a > exp_score_b or
                    (
                        exp_score_a == exp_score_b and
                        expansion_len_a > expansion_len_b
                    )
                ):
                    slice_length_a += expansion_len_a
                    slice_length_b += len(use_b_string[:slice_start_b])
                    slice_start_a -= expansion_len_a
                    slice_start_b = 0
                else:
                    slice_length_a += len(bgc_a.gene_string[:slice_start_a])
                    slice_length_b += expansion_len_b
                    slice_start_a = 0
                    slice_start_b -= expansion_len_b

            else:
                # A is shorter upstream. Assume complete extension. Find B's
                # extension
                if slice_start_a < slice_start_b:
                    exp_score_b, exp_length_b = score_expansion(
                        use_b_string[:slice_start_b],
                        bgc_a.gene_string[:slice_start_a],
                        False
                    )

                    slice_length_a += len(bgc_a.gene_string[:slice_start_a])
                    slice_length_b += exp_length_b
                    slice_start_a = 0
                    slice_start_b -= exp_length_b
                else:
                    exp_score_a, exp_length_a = score_expansion(
                        bgc_a.gene_string[:slice_start_a],
                        use_b_string[:slice_start_b],
                        False
                    )

                    slice_length_a += exp_length_a
                    slice_length_b += len(use_b_string[:slice_start_b])
                    slice_start_a -= exp_length_a
                    slice_start_b = 0

            # RIGHT SIDE
            # check which side is the shortest downstream. If both BGCs have
            # the same length left, choose the one with the best expansion
            downstream_a = bgc_a.num_genes - slice_start_a - slice_length_a
            downstream_b = bgc_b.num_genes - slice_start_b - slice_length_b
            if downstream_a == downstream_b:
                # assume complete extension of A, try to expand B
                exp_score_b, exp_len_b = score_expansion(
                    use_b_string[slice_start_b+slice_length_b:],
                    bgc_a.gene_string[slice_start_a+slice_length_a:],
                    True
                )
                # assume complete extension of B, try to expand A
                exp_score_a, exp_len_a = score_expansion(
                    bgc_a.gene_string[slice_start_a+slice_length_a:],
                    use_b_string[slice_start_b+slice_length_b:],
                    True
                )


                if (
                    (exp_score_a == exp_score_b and exp_len_a > exp_len_b) or
                    exp_score_a > exp_score_b
                ):
                    slice_length_a += exp_len_a
                    b_substring = use_b_string[slice_start_b+slice_length_b:]
                    slice_length_b += len(b_substring)
                else:
                    a_substring = bgc_a.gene_string[slice_start_a+slice_length_a:]
                    slice_length_a += len(a_substring)
                    slice_length_b += exp_len_b

            else:
                if downstream_a < downstream_b:
                    # extend all of remaining A
                    exp_score_b, exp_len_b = score_expansion(
                        use_b_string[slice_start_b+slice_length_b:],
                        bgc_a.gene_string[slice_start_a+slice_length_a:],
                        True
                    )

                    a_substring = bgc_a.gene_string[slice_start_a+slice_length_a:]
                    slice_length_a += len(a_substring)
                    slice_length_b += exp_len_b

                else:
                    exp_score_a, exp_len_a = score_expansion(
                        bgc_a.gene_string[slice_start_a+slice_length_a:],
                        use_b_string[slice_start_b+slice_length_b:],
                        True
                    )
                    slice_length_a += exp_len_a
                    slice_length_b += len(use_b_string[slice_start_b+slice_length_b:])

        if min(slice_length_a, slice_length_b) < 5:
            return

        # First test passed. Find if there is a biosynthetic gene in both slices
        # (note that even if they are, currently we don't check whether it's
        # actually the _same_ gene)
        biosynthetic_hit_a = False
        biosynthetic_hit_b = False

        for biosynthetic_position in bgc_a.bio_synth_core_positions:
            if (
                biosynthetic_position >= slice_start_a and
                biosynthetic_position <= slice_start_a + slice_length_a
            ):
                biosynthetic_hit_a = True
                break


        # return to original orientation if needed
        if reverse:
            slice_start_b = bgc_b.num_genes - slice_start_b - slice_length_b

        # using original bio_synth_core_positions
        for biosynthetic_position in bgc_b.bio_synth_core_positions:
            if biosynthetic_position >= slice_start_b and biosynthetic_position <= (slice_start_b + slice_length_b):
                biosynthetic_hit_b = True
                break

        # finally...
        if biosynthetic_hit_a and biosynthetic_hit_b:
            self.a_dom_start = sum(bgc_a.gene_domain_counts[:slice_start_a])
            self.a_dom_end = self.a_dom_start + sum(bgc_a.gene_domain_counts[slice_start_a:slice_start_a+slice_length_a])
            self.a_dom_set = set(bgc_a.ordered_domain_list[self.a_dom_start:self.a_dom_end])


            self.b_dom_start = sum(bgc_b.gene_domain_counts[:slice_start_b])
            self.b_dom_end = self.b_dom_start + sum(bgc_b.gene_domain_counts[slice_start_b:slice_start_b+slice_length_b])
            self.b_dom_set = set(bgc_b.ordered_domain_list[self.b_dom_start:self.b_dom_end])

            self.intersect = self.a_dom_set & self.b_dom_set

            # re-adjust the indices for each domain so we get only the sequence
            # tags in the selected slice. First step: find out which is the
            # first copy of each domain we're using
            for domain in bgc_a.ordered_domain_list[:self.a_dom_start]:
                self.a_dom_seq_slice_bot[domain] += 1
            for domain in bgc_b.ordered_domain_list[:self.b_dom_start]:
                self.b_dom_seq_slice_bot[domain] += 1

            # Step 2: work with the last copy of each domain.
            # Step 2a: make top = bottom
            for domain in self.a_dom_set:
                self.a_dom_seq_slice_top[domain] = self.a_dom_seq_slice_bot[domain]
            for domain in self.b_dom_set:
                self.b_dom_seq_slice_top[domain] = self.b_dom_seq_slice_bot[domain]

            # Step 2b: increase top with the domains in the slice
            for domain in bgc_a.ordered_domain_list[self.a_dom_start:self.a_dom_end]:
                self.a_dom_seq_slice_top[domain] += 1
            for domain in bgc_b.ordered_domain_list[self.b_dom_start:self.b_dom_end]:
                self.b_dom_seq_slice_top[domain] += 1
