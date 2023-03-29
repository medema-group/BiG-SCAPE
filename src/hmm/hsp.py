"""Contains class to describe domain hits and alignments"""

# from python
from __future__ import annotations

# from dependencies
from pyhmmer.plan7 import Domain

# from other modules
from src.genbank import CDS


class HSP:
    """Describes a CDS - Domain relationship"""

    def __init__(self, cds: CDS, domain: Domain, score) -> None:
        self.cds = cds
        self.domain = domain
        self.score = score

    def save(self):
        """Saves this object to a database"""
        pass

    def has_overlap(self, hsp_b: HSP):
        """Return True if there is overlap between two regions"""
        # a left of b
        if self.cds.nt_stop < hsp_b.cds.nt_start:
            return False
        # b left of a
        if hsp_b.cds.nt_stop < self.cds.nt_start:
            return False

        # all other cases should have overlap
        return True

    def len_overlap(self, hsp_b: HSP):
        """Returns the length of an overlapping sequence"""

        if self.cds.nt_start < hsp_b.cds.nt_start:
            left = hsp_b.cds.nt_start
        else:
            left = self.cds.nt_start

        if self.cds.nt_stop > hsp_b.cds.nt_stop:
            right = hsp_b.cds.nt_stop
        else:
            right = self.cds.nt_stop

        # limit to > 0
        return max(0, right - left)

    def __repr__(self) -> str:
        return ",".join(
            map(
                str,
                [
                    self.cds.nt_start,
                    self.cds.nt_stop,
                    self.cds.strand,
                    self.domain,
                    self.score,
                ],
            )
        )

    def __eq__(self, __o: object) -> bool:
        """Return whether this HSP object and another object are equal

        Will always return false if the object compared is not an instance of HSP

        Args:
            __o (object): Target object to compare to

        Returns:
            bool: True if __o and self are equal
        """
        if not isinstance(__o, HSP):
            return False

        return all(
            [
                __o.cds.nt_start == self.cds.nt_start,
                __o.cds.nt_stop == self.cds.nt_stop,
                __o.domain == self.domain,
                __o.score == self.score,
            ]
        )


def hsp_overlap_filter(hsp_list: list[HSP], overlap_cutoff=0.1) -> list[HSP]:
    """Filters overlapping HSP

    Args:
        hsp_list (list[HSP]): a list of high scoring protein hits
        overlap_cutoff: The maximum percentage sequence length domains can overlap
            without being discarded

    Returns:
        list[HSP]: a filtered list of high scoring protein hits
    """
    # for this the hsps have to be sorted by any of the start coordinates
    # this is for the rare case where two hsps have the same bit score.
    # to replicate BiG-SCAPE output, the hsp with a lower start coordinate
    # will end up being chosen
    sorted_hsps = sorted(hsp_list, key=lambda hsp: hsp.cds.nt_start)

    new_list = []
    for i in range(len(sorted_hsps) - 1):
        for j in range(i + 1, len(sorted_hsps)):
            hsp_a: HSP = sorted_hsps[i]
            hsp_b: HSP = sorted_hsps[j]

            # check if we are the same CDS
            if hsp_a == hsp_b:
                # check if there is overlap between the domains
                if not hsp_a.has_overlap(hsp_b):
                    continue

                len_aa_overlap = hsp_a.len_overlap(hsp_b)
                overlap_perc_loc1 = len_aa_overlap / (
                    hsp_a.cds.nt_stop - hsp_a.cds.nt_start
                )
                overlap_perc_loc2 = len_aa_overlap / (
                    hsp_b.cds.nt_stop - hsp_b.cds.nt_start
                )
                # check if the amount of overlap is significant
                if (
                    overlap_perc_loc1 <= overlap_cutoff
                    and overlap_perc_loc2 <= overlap_cutoff
                ):
                    continue
                # rounding with 1 decimal to mirror domtable files
                hsp_a_score = round(float(hsp_a.score), 1)
                hsp_b_score = round(float(hsp_b.score), 1)
                if hsp_a_score >= hsp_b_score:  # see which has a better score
                    new_list.append(hsp_a)
                elif hsp_a_score < hsp_b_score:
                    new_list.append(hsp_b)

    return new_list


class HSPAlignment:
    """Describes a HSP - Domain alignment"""

    def __init__(self, hsp: HSP, domain: str, alignment: str) -> None:
        self.hsp = hsp
        self.domain = domain
        self.alignment = alignment

    def save(self):
        """Saves this object to a database"""
        pass

    def __repr__(self) -> str:
        return f"Alignment of hsp {str(self.hsp)} and domain {self.domain}: {self.alignment}"

    def __eq__(self, __o: object) -> bool:
        """Return whether this HSP alignment and another object are equal

        Will always return false if the object compared is not an instance of
        HSPAlignment

        Args:
            __o (object): Target object to compare to

        Returns:
            bool: True if __o and self are equal
        """
        if not isinstance(__o, HSPAlignment):
            return False

        conditions = [
            self.hsp == __o.hsp,
            self.domain == __o.domain,
            self.alignment == __o.alignment,
        ]
        return all(conditions)
