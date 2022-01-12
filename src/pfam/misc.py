import sys
import os


def get_domain_list(filename):
    """Convert the Pfam string in the .pfs files to a Pfam list"""

    domains = []
    with open(filename, "r") as handle:
        # note: cannot filter the domain version as the BGCs dictionary needs the
        # complete Pfam id
        domains = handle.readline().strip().split(" ")

    return domains



def check_overlap(pfd_matrix, overlap_cutoff):
    """Check if domains overlap for a certain overlap_cutoff.
     If so, remove the domain(s) with the lower score."""

    delete_list = []
    for i in range(len(pfd_matrix)-1):
        for j in range(i+1, len(pfd_matrix)):
            row1 = pfd_matrix[i]
            row2 = pfd_matrix[j]

            #check if we are the same CDS
            if row1[-1] == row2[-1]:
                #check if there is overlap between the domains
                if not no_overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4])):
                    overlapping_aminoacids = len_overlap(int(row1[3]), int(row1[4]), int(row2[3]), int(row2[4]))
                    overlap_perc_loc1 = overlap_perc(overlapping_aminoacids, int(row1[4])-int(row1[3]))
                    overlap_perc_loc2 = overlap_perc(overlapping_aminoacids, int(row2[4])-int(row2[3]))
                    #check if the amount of overlap is significant
                    if overlap_perc_loc1 > overlap_cutoff or overlap_perc_loc2 > overlap_cutoff:
                        if float(row1[1]) >= float(row2[1]): #see which has a better score
                            delete_list.append(row2)
                        elif float(row1[1]) < float(row2[1]):
                            delete_list.append(row1)

    for lst in delete_list:
        try:
            pfd_matrix.remove(lst)
        except ValueError:
            pass

    # for some reason, some coordinates in genbank files have ambiguous
    # starting/ending positions for the CDS. In this case we need the
    # absolute start position of each domain so the loci start coordinate
    # needs to be checked
    absolute_start_positions = [] # we have to take care of strandedness of the origin gene
    for row in pfd_matrix:
        if row[7][0] == "<" or row[7][0] == ">":
            row[7] = row[7][1:]
        if row[8][0] == "<" or row[8][0] == ">":
            row[8] = row[8][1:]

        loci_start = int(row[7])
        loci_end = int(row[8])
        domain_start = int(row[3])
        domain_end = int(row[4])

        # uses nucleotide coordinates for sorting
        width = 3*(domain_end - domain_start)
        strand = row[9].split(":")[-1]
        if strand == "+":
            absolute_start_positions.append(3*domain_start + loci_start)
        else:
            absolute_start_positions.append(loci_end - 3*domain_start - width)

    try:
        # sorted by absolute position: env coordinate start + gene locus start
        pfd_matrix = [x for (y, x) in sorted(zip(absolute_start_positions, pfd_matrix), key=lambda pair: pair[0])]
    except ValueError as e:
        print("Something went wrong during the sorting of the fourth column (ValueError): " + str(e))
        print(pfd_matrix)
        sys.exit()

    domains = []
    for row in pfd_matrix:
        # removing the domain version at this point is not a very good idea
        # because we need it complete for hmmalign
        domains.append(row[5]) #save the pfam domains for the .pfs file

    return pfd_matrix, domains


def no_overlap(a_start, a_end, b_start, b_end):
    """Return True if there is no overlap between two regions"""
    if a_start < b_start and a_end < b_start:
        return True
    elif a_start > b_end and a_end > b_end:
        return True
    else:
        return False


def overlap_perc(overlap, len_seq):
    return float(overlap) / len_seq


def len_overlap(a_start, a_end, b_start, b_end):
    """Returns the amount of overlapping nucleotides"""

    if a_start < b_start:
        cor1 = a_start
    else:
        cor1 = b_start

    if a_end > b_end:
        cor2 = a_end
    else:
        cor2 = b_end

    total_region = cor2 - cor1
    sum_len = (a_end - a_start) + (b_end - b_start)

    return sum_len - total_region

def generate_pfam_colors_matrix(pfam_domain_colors):
    '''

    :param pfam_domain_colors: tab-delimited file
    :return: dictionary with pfam ID as key and rgb colors as value
    '''
    pfam_colors = {}

    if os.path.isfile(pfam_domain_colors):
        print("  Found file with Pfam domain colors")
        with open(pfam_domain_colors, "r") as cat_handle:
            for line in cat_handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    row = line.strip().split("\t")
                    domain = row[0]
                    rgb = row[-1]
                    pfam_colors[domain] = rgb
    else:
        print("  File pfam_domain_colors was NOT found")

    return pfam_colors
