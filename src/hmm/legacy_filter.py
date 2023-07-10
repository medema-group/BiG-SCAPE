"""Contains legacy (BiG-SCAPE v1.0) code to filter and remove HSP hits depending on the
overlap with other hits
"""


# from other modules
from src.hmm import HSP


def legacy_filter_overlap(hsp_list: list[HSP], overlap_cutoff) -> list[HSP]:
    """Perform filtering of HSPs based on the 1.0 approach of all-vs-all comparison

    Args:
        hsp_list (list[HSP]): List of HSPs to filter
        overlap_cutoff (float): cutoff percentage of total overlapping nucleotides
    Returns:
        (list[HSP]): A new, filtered list of HSPs
    """

    delete_list = []
    new_list = hsp_list.copy()

    for i in range(len(hsp_list) - 1):
        for j in range(i + 1, len(hsp_list)):
            row1 = hsp_list[i]
            row2 = hsp_list[j]

            # check if we are the same CDS
            if row1.cds == row2.cds:
                # check if there is overlap between the domains
                if not no_overlap(
                    row1.env_start, row1.env_stop, row2.env_start, row2.env_stop
                ):
                    overlapping_aminoacids = overlap(
                        int(row1.env_start),
                        int(row1.env_stop),
                        int(row2.env_start),
                        int(row2.env_stop),
                    )
                    overlap_perc_loc1 = overlap_perc(
                        overlapping_aminoacids, int(row1.env_stop) - int(row1.env_start)
                    )
                    overlap_perc_loc2 = overlap_perc(
                        overlapping_aminoacids, int(row2.env_stop) - int(row2.env_start)
                    )
                    # check if the amount of overlap is significant
                    if (
                        overlap_perc_loc1 > overlap_cutoff
                        or overlap_perc_loc2 > overlap_cutoff
                    ):
                        if float(row1.score) >= float(
                            row2.score
                        ):  # see which has a better score
                            delete_list.append(row2)
                        elif float(row1.score) < float(row2.score):
                            delete_list.append(row1)

    for lst in delete_list:
        try:
            new_list.remove(lst)
        except ValueError:
            pass
    return new_list


def no_overlap(locA1, locA2, locB1, locB2):
    """Return True if there is no overlap between two regions"""
    if locA1 < locB1 and locA2 < locB1:
        return True
    elif locA1 > locB2 and locA2 > locB2:
        return True
    else:
        return False


def overlap_perc(overlap, len_seq):
    return float(overlap) / len_seq


def overlap(a_start, a_stop, b_start, b_stop):
    """Returns the amount of overlapping nucleotides"""

    if a_start < b_start:
        overlap_start = a_start
    else:
        overlap_start = b_start

    if a_stop > b_stop:
        overlap_stop = a_stop
    else:
        overlap_stop = b_stop

    total_region = overlap_stop - overlap_start
    sum_len = (a_stop - a_start) + (b_stop - b_start)

    return sum_len - total_region
