"""Contains methods to calculate the comparable region of a BGC pair, which is used to
determine the region between two BGCs for which to calculate AI and DSS
"""

# from python
from difflib import SequenceMatcher

# from this module
from .binning import BGCPair


def find_dom_list_orientation():
    raise NotImplementedError()


# TODO: is this all we need as output? or do we really still need all the start and stop
# stuff?
def get_pair_domain_lcs(pair: BGCPair) -> tuple[list[str], bool]:
    """Retrieve the longest common subsequence of domains for a pair of BGC records
    Also returns whether the sequence is reversed

    Args:
        pair (BGCPair): Pair of BGCs to retrieve the LCS for

    Returns:
        tuple[list[str], bool]: List of strings that represents the LCS of domains and
        a boolean indicating whether the LCS is was found in the reverse sequence
    """
    a_hsps = pair.region_a.get_hsps()
    b_hsps = pair.region_b.get_hsps()

    a_domains = [a_hsp.domain for a_hsp in a_hsps]
    b_domains = [b_hsp.domain for b_hsp in b_hsps]

    # forward
    seqmatch = SequenceMatcher(None, a_domains, b_domains)
    match = seqmatch.find_longest_match(0, len(a_domains), 0, len(b_domains))
    a_start_fwd = match[0]
    fwd_match_len = match[2]

    # reverse
    seqmatch = SequenceMatcher(None, a_domains, b_domains[::-1])
    match = seqmatch.find_longest_match(0, len(a_domains), 0, len(b_domains))
    a_start_rev = match[0]
    rev_match_len = match[2]

    if fwd_match_len >= rev_match_len:
        reverse = False
        a_start = a_start_fwd
        match_len = fwd_match_len
    else:
        reverse = True
        a_start = a_start_rev
        match_len = rev_match_len

    if match_len == 0:
        return ([], reverse)

    a_stop = a_start + match_len

    return a_domains[a_start:a_stop], reverse
