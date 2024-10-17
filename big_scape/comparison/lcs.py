"""Contains functions for computing the longest common subsequence between two lists
of domains or CDS in a RecordPair

Into any of the main functions in this module goes a RecordPair object, which contains
two regions. These regions can be either protoclusters or full regions. The functions
in this module are used to find the longest common subsequence between the two regions
in the RecordPair. This is used in extension as a "seed" for the extension.

The result from the main functions are tuples with the following structure:
    (a_start, a_stop, b_start, b_stop, reverse)
Where a and b are the two regions in the RecordPair. Start is inclusive, stop is
exclusive. Reverse is a boolean indicating whether the match is in reverse

NOTE: The matches correspond to slices of region CDS that do not have domains!
"""

# from python
import logging
from typing import Any
from difflib import Match, SequenceMatcher
from sqlalchemy import select, update, or_, and_

# from other modules
import big_scape.genbank as bs_genbank
import big_scape.comparison.record_pair as bs_comparison
import big_scape.hmm as bs_hmm
from big_scape.data import DB


def find_lcs(list_a: list[Any], list_b: list[Any]) -> tuple[Match, list[Match]]:
    """Detect longest common substring using sequencematcher

    Args:
        list_a (list[T]): A list of hashable objects
        list_b (list[T]): Second list of hashable objects

    Returns:
        tuple[int, int, int]: start, stop and length of longest common substring
    """
    seqmatch = SequenceMatcher(None, list_a, list_b, False)
    match = seqmatch.find_longest_match(0, len(list_a), 0, len(list_b))
    matching_blocks = seqmatch.get_matching_blocks()
    return match, matching_blocks


def find_bio_or_middle_lcs(
    a_cds: list[bs_genbank.CDS],
    b_cds: list[bs_genbank.CDS],
    a_domains: list[bs_hmm.HSP],
    b_domains: list[bs_hmm.HSP],
    matching_blocks: list[Match],
    matching_blocks_rev: list[Match],
    a_domain_cds_idx: dict[int, int],
    b_domain_cds_idx: dict[int, int],
) -> tuple[int, int, int, int, bool]:
    """Find the most central match out of all LCS matches, or a match containing a
    biosynthetic gene

    This is done by first selecting the shorter record in terms of CDS, and then
    finding the match that is closest to the middle of the CDS

    Args:
        a_cds (list[CDS]): List of CDS for A
        b_cds (list[CDS]): List of CDS for B
        a_domains (list[HSP]): List of domains for A
        b_domains (list[HSP]): List of domains for B
        matching_blocks (list[Match]): List of matching blocks
        matching_blocks_rev (list[Match]): List of matching blocks in reverse
        a_domain_cds_idx (dict[int, int]): Dictionary of domain idx to cds idx for A
        b_domain_cds_idx (dict[int, int]): Dictionary of domain idx to cds idx for B

    Returns:
        tuple[int, int, int, int, bool]: a_start, a_stop, b_start, b_stop, reverse
    """

    # first find which region is shorter in terms of cds
    a_cds_len = len(a_cds)
    b_cds_len = len(b_cds)

    # default to A. if B is shorter, use B
    if a_cds_len <= b_cds_len:
        use_cds = a_cds
        use_domains = a_domains
        matching_block_idx = 0
    else:
        use_cds = b_cds
        use_domains = b_domains
        matching_block_idx = 1

    # generate a CDS to index dict
    cds_idx_dict = {cds: i for i, cds in enumerate(use_cds)}

    # go through all LCS matches and find the one with the most central CDS
    middle = len(use_cds) / 2
    min = None
    for matching_block in matching_blocks:
        # I don't even know why there is a match of len 0 when there are matches
        # of len 1
        if matching_block[2] == 0:
            continue

        idx = matching_block[matching_block_idx]

        domain = use_domains[idx]
        cds_idx = cds_idx_dict[domain.cds]

        # find the distance to the middle
        distance = abs(middle - cds_idx)

        if min is None or distance < min:
            min = distance
            a_start = matching_block[0]
            a_stop = matching_block[0] + matching_block[2]
            b_start = matching_block[1]
            b_stop = matching_block[1] + matching_block[2]
            reverse = False

    # go through all reverse, too
    for matching_block in matching_blocks_rev:
        if matching_block[2] == 0:
            continue

        idx = matching_block[matching_block_idx]

        domain = use_domains[idx]
        cds_idx = cds_idx_dict[domain.cds]

        # find the distance to the middle
        distance = abs(middle - cds_idx)

        if min is None or distance < min:
            min = distance
            a_start = matching_block[0]
            a_stop = matching_block[0] + matching_block[2]
            b_start = len(b_domains) - matching_block[1] - matching_block[2]
            b_stop = len(b_domains) - matching_block[1]
            reverse = True

    a_cds_start = a_domain_cds_idx[a_start]
    a_cds_stop = a_domain_cds_idx[a_stop - 1] + 1
    b_cds_start = b_domain_cds_idx[b_start]
    b_cds_stop = b_domain_cds_idx[b_stop - 1] + 1

    return a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse


def find_domain_lcs_region(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two lists of domains

    This takes CDS as arguments, but uses the domains within the CDS to find the LCS

    NOTE: The LCS correspond to slices of region CDS that have domains only!
    These slices need to be converted to full CDS ranges later

    Approach:
    - Pick the largest LCS with a biosynthetic gene
    - If there are multiple LCS with a biosynthetic gene of equal length, pick the one
        that is closest to the middle of the region
    - If there are no LCS with a biosynthetic gene, pick the largest
    - If there are multiple LCS of equal length, pick the one that is closest to the
        middle of the region

    Args:
        a_cds (list[CDS]): List of CDS
        b_cds (list[CDS]): List of CDS

    Returns:
        tuple[int, int, int, int, int, int, int, int, bool]: a_start, a_stop,
            b_start, b_stop, a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse
    """
    logging.debug("region lcs")

    # these are regions, so we can get the full range of CDS
    a_cds = list(pair.record_a.get_cds_with_domains())
    b_cds = list(pair.record_b.get_cds_with_domains())

    # working on domains, not cds
    a_domains: list[bs_hmm.HSP] = []
    b_domains: list[bs_hmm.HSP] = []

    # dictionary of domain index to cds index to quickly find the cds of a domain
    a_domain_cds_idx = {}
    b_domain_cds_idx = {}

    # list of domain idx whose genes are biosynthetic, to quickly find if a domain is
    # part of a biosynthetic gene
    a_biosynthetic_domain_cds: list[int] = []
    b_biosynthetic_domain_cds: list[int] = []

    for cds_idx, cds in enumerate(a_cds):
        for i in range(len(a_domains), len(a_domains) + len(cds.hsps)):
            a_domain_cds_idx[i] = cds_idx

        if cds.strand == 1:
            a_domains.extend(cds.hsps)
        elif cds.strand == -1:
            a_domains.extend(cds.hsps[::-1])

        if cds.gene_kind == "biosynthetic":
            a_biosynthetic_domain_cds.append(cds_idx)

    for cds_idx, cds in enumerate(b_cds):
        for i in range(len(b_domains), len(b_domains) + len(cds.hsps)):
            b_domain_cds_idx[i] = cds_idx

        if cds.strand == 1:
            b_domains.extend(cds.hsps)
        elif cds.strand == -1:
            b_domains.extend(cds.hsps[::-1])

        if cds.gene_kind == "biosynthetic":
            b_biosynthetic_domain_cds.append(cds_idx)

    # forward
    match, matching_blocks_fwd = find_lcs(a_domains, b_domains)
    fwd_match_len = match[2]

    # reverse
    match, matching_blocks_rev = find_lcs(a_domains, b_domains[::-1])
    rev_match_len = match[2]

    # quickly check if we didn't find an LCS
    if fwd_match_len == 0 and rev_match_len == 0:
        logging.error(
            "No match found in LCS. This should not happen after first jaccard"
        )
        logging.error("a domains: %s", a_domains)
        logging.error("b domains: %s", b_domains)
        raise RuntimeError("No match found in LCS.")

    # now we need to do something silly. we want to assemble a list of these matching
    # blocks, but we want to keep track of whether they are in reverse or not. this will
    # make it so we can do everything we need to do in one loop later
    matching_block_dirs = []
    for matching_block in matching_blocks_fwd:
        matching_block_dirs.append((matching_block + (False,)))

    for matching_block in matching_blocks_rev:
        matching_block_dirs.append((matching_block + (True,)))

    # this is where the fun begins. we will use these lists to decide which match to
    # return later

    # tuple is idx, length, reverse
    longest_biosynthetic: list[tuple[int, int, bool]] = []
    longest: list[tuple[int, int, bool]] = []
    # tuple is idx, distance to middle, reverse
    central_biosynthetic: list[tuple[int, int, bool]] = []
    central: list[tuple[int, int, bool]] = []

    for match_idx, matching_block_dir in enumerate(matching_block_dirs):
        start_a = matching_block_dir[0]
        stop_a = matching_block_dir[0] + matching_block_dir[2]
        start_b = matching_block_dir[1]
        stop_b = matching_block_dir[1] + matching_block_dir[2]
        length = matching_block_dir[2]
        reverse = matching_block_dir[3]

        # I don't understand why, but zero-length blocks exist sometimes. skip them
        if length == 0:
            continue

        # check if the match contains a biosynthetic gene
        has_biosynthetic = False
        for biosynthetic_idx in a_biosynthetic_domain_cds:
            if (
                a_domain_cds_idx[start_a]
                <= biosynthetic_idx
                < a_domain_cds_idx[stop_a - 1] + 1
            ):
                has_biosynthetic = True
                break

        for biosynthetic_idx in b_biosynthetic_domain_cds:
            if (
                b_domain_cds_idx[start_b]
                <= biosynthetic_idx
                < b_domain_cds_idx[stop_b - 1] + 1
            ):
                has_biosynthetic = True
                break

        # select for biosynthetic or normal list
        use_longest_list = longest_biosynthetic if has_biosynthetic else longest
        use_central_list = central_biosynthetic if has_biosynthetic else central

        # Length

        # clear the list if it's not empty and the current match is longer
        if len(use_longest_list) > 0 and length > use_longest_list[0][1]:
            use_longest_list.clear()

        # add the match to the list if it's empty or the current match is equal length
        # or longer than the existing match
        if len(use_longest_list) == 0 or length >= use_longest_list[0][1]:
            use_longest_list.append((match_idx, length, reverse))

        # distance to middle

        # use the shorter cds list to determine the distance to middle
        use_cds = a_cds if len(a_cds) <= len(b_cds) else b_cds
        use_idx = a_domain_cds_idx if len(a_cds) <= len(b_cds) else b_domain_cds_idx
        use_start = start_a if len(a_cds) <= len(b_cds) else start_b
        use_stop = stop_a if len(a_cds) <= len(b_cds) else stop_b

        middle = len(use_cds) / 2
        # calculate the distance from either side of the match to the middle
        distance = round(
            min(
                abs(middle - use_idx[use_start]),
                abs(middle - use_idx[use_stop - 1]),
            )
        )

        # clear the list if it's not empty and the current match is closer to the middle
        if len(use_central_list) > 0 and distance < use_central_list[0][1]:
            use_central_list.clear()

        # add the match to the list if it's empty or the current match is equal distance
        # or closer to the middle than the existing match
        if len(use_central_list) == 0 or distance <= use_central_list[0][1]:
            use_central_list.append((match_idx, distance, reverse))

    # now we have everything we need. we need to decide which match to return
    # just go top to bottom and return the first match in the list
    # remember that the lists are [match_idx, length/dist, reverse]
    # refer to docstring for decision making here
    if len(longest_biosynthetic) == 1:
        match_idx = longest_biosynthetic[0][0]

    elif len(central_biosynthetic) > 0:
        match_idx = central_biosynthetic[0][0]

    elif len(longest) == 1:
        match_idx = longest[0][0]

    elif len(central) > 0:
        match_idx = central[0][0]

    else:
        # this should never happen
        raise RuntimeError("No match found in LCS.")

    # match bounds on domain level will be inclusive to ensure a consistent exclusive
    # stop on cds level
    relevant_match = matching_block_dirs[match_idx]
    a_start = relevant_match[0]
    a_stop = relevant_match[0] + relevant_match[2] - 1  # inclusive stop
    b_start = relevant_match[1]
    b_stop = relevant_match[1] + relevant_match[2] - 1  # inclusive stop
    reverse = relevant_match[3]

    a_cds_start = a_domain_cds_idx[a_start]
    # cds stop may be end of cds
    if a_stop == len(a_domains) - 1:
        a_cds_stop = len(a_cds)
    else:
        a_cds_stop = a_domain_cds_idx[a_stop] + 1  # exclusive stop

    # fix b start and stop if in reverse. this means first flipping the domain indexes, getting the
    # cds index, and then flipping the cds index. yay.
    if reverse:
        old_start = b_start
        old_stop = b_stop
        b_start = len(b_domains) - b_stop - 1
        b_stop = len(b_domains) - old_start - 1  # inclusive stop

    b_cds_start = b_domain_cds_idx[b_start]

    # cds stop may be end of cds
    if b_stop == len(b_domains) - 1:
        b_cds_stop = len(b_cds)
    else:
        b_cds_stop = b_domain_cds_idx[b_stop] + 1  # exclusive stop

    if reverse:
        old_cds_start = b_cds_start
        b_cds_start = len(b_cds) - b_cds_stop
        b_cds_stop = len(b_cds) - old_cds_start  # exclusive stop

    # final check: it could happen that the start and stop of the domain LCS is in the
    # same CDS. in this case, the stop needs to be incremented by one
    if a_cds_start == a_cds_stop:
        a_cds_stop += 1

    if b_cds_start == b_cds_stop:
        b_cds_stop += 1

    # revert inclusive domain stop and index flipping before returning
    if reverse:
        b_start = old_start
        b_stop = old_stop
    a_stop += 1
    b_stop += 1

    return (
        a_start,
        a_stop,
        b_start,
        b_stop,
        a_cds_start,
        a_cds_stop,
        b_cds_start,
        b_cds_stop,
        reverse,
    )


def find_domain_lcs_protocluster(
    pair: bs_comparison.RecordPair,
) -> tuple[int, int, int, int, int, int, int, int, bool]:
    """Find the longest stretch of matching domains between two protocluster records,
    using domains

    NOTE: The LCS correspond to slices of region CDS that have domains only!
    These slices need to be converted to full CDS ranges later

    Args:
        pair (RecordPair): RecordPair object

    Returns:
        tuple[int, int, int, int, int, int, int, int, bool]: a_start, a_stop,
            b_start, b_stop, a_cds_start, a_cds_stop, b_cds_start, b_cds_stop, reverse
    """
    logging.debug("pc lcs")

    # we really need protoclusters here
    if not isinstance(pair.record_a, bs_genbank.ProtoCluster):
        raise TypeError("record_a must be a protocluster")

    if not isinstance(pair.record_b, bs_genbank.ProtoCluster):
        raise TypeError("record_b must be a protocluster")

    a_cds = list(pair.record_a.get_cds_with_domains())
    b_cds = list(pair.record_b.get_cds_with_domains())

    # working on domains, not cds
    a_domains: list[bs_hmm.HSP] = []
    b_domains: list[bs_hmm.HSP] = []

    # dictionary of domain index to cds index to quickly find the cds of a domain
    a_domain_cds_idx = {}
    b_domain_cds_idx = {}

    for cds_idx, cds in enumerate(a_cds):
        for i in range(len(a_domains), len(a_domains) + len(cds.hsps)):
            a_domain_cds_idx[i] = cds_idx

        if cds.strand == 1:
            a_domains.extend(cds.hsps)
        elif cds.strand == -1:
            a_domains.extend(cds.hsps[::-1])

    for cds_idx, cds in enumerate(b_cds):
        for i in range(len(b_domains), len(b_domains) + len(cds.hsps)):
            b_domain_cds_idx[i] = cds_idx

        if cds.strand == 1:
            b_domains.extend(cds.hsps)
        elif cds.strand == -1:
            b_domains.extend(cds.hsps[::-1])

    # forward
    match, matching_blocks_fwd = find_lcs(a_domains, b_domains)
    fwd_match_len = match[2]

    # reverse
    match, matching_blocks_rev = find_lcs(a_domains, b_domains[::-1])
    rev_match_len = match[2]

    # quickly check if we didn't find an LCS
    if fwd_match_len == 0 and rev_match_len == 0:
        logging.error(
            "No match found in LCS. This should not happen after first jaccard"
        )
        logging.error("a domains: %s", a_domains)
        logging.error("b domains: %s", b_domains)
        raise RuntimeError("No match found in LCS.")

    # now we need to do something silly. we want to assemble a list of these matching
    # blocks, but we want to keep track of whether they are in reverse or not. this will
    # make it so we can do everything we need to do in one loop later
    matching_block_dirs = []
    for matching_block in matching_blocks_fwd:
        matching_block_dirs.append((matching_block + (False,)))

    for matching_block in matching_blocks_rev:
        matching_block_dirs.append((matching_block + (True,)))

    # tuple is idx, length, reverse
    longest_protocore: list[tuple[int, int, bool]] = []
    longest: list[tuple[int, int, bool]] = []
    # tuple is idx, distance to middle, reverse
    central_protocore: list[tuple[int, int, bool]] = []
    central: list[tuple[int, int, bool]] = []

    for match_idx, matching_block_dir in enumerate(matching_block_dirs):
        start_a = matching_block_dir[0]
        stop_a = matching_block_dir[0] + matching_block_dir[2]
        start_b = matching_block_dir[1]
        stop_b = matching_block_dir[1] + matching_block_dir[2]
        length = matching_block_dir[2]
        reverse = matching_block_dir[3]

        # I don't understand why, but zero-length blocks exist sometimes. skip them
        if length == 0:
            continue

        # check if this match is (partly) in the protocore
        in_protocore = False
        for protocore_idx in pair.record_a.proto_core_cds_idx:
            if start_a <= protocore_idx < stop_a:
                in_protocore = True
                break

        for protocore_idx in pair.record_b.proto_core_cds_idx:
            if start_b <= protocore_idx < stop_b:
                in_protocore = True
                break

        # select for biosynthetic or normal list
        use_longest_list = longest_protocore if in_protocore else longest
        use_central_list = central_protocore if in_protocore else central

        # Length

        # clear the list if it's not empty and the current match is longer
        if len(use_longest_list) > 0 and length > use_longest_list[0][1]:
            use_longest_list.clear()

        # add the match to the list if it's empty or the current match is equal length
        # or longer than the existing match
        if len(use_longest_list) == 0 or length >= use_longest_list[0][1]:
            use_longest_list.append((match_idx, length, reverse))

        # distance to middle

        # use the shorter cds list to determine the distance to middle
        use_cds = a_cds if len(a_cds) <= len(b_cds) else b_cds
        use_idx = a_domain_cds_idx if len(a_cds) <= len(b_cds) else b_domain_cds_idx
        use_start = start_a if len(a_cds) <= len(b_cds) else start_b
        use_stop = stop_a if len(a_cds) <= len(b_cds) else stop_b

        middle = len(use_cds) / 2
        # calculate the distance from either side of the match to the middle
        distance = round(
            min(
                abs(middle - use_idx[use_start]),
                abs(middle - use_idx[use_stop - 1]),
            )
        )

        # clear the list if it's not empty and the current match is closer to the middle
        if len(use_central_list) > 0 and distance < use_central_list[0][1]:
            use_central_list.clear()

        # add the match to the list if it's empty or the current match is equal distance
        # or closer to the middle than the existing match
        if len(use_central_list) == 0 or distance <= use_central_list[0][1]:
            use_central_list.append((match_idx, distance, reverse))

    # now we have everything we need. we need to decide which match to return
    # just go top to bottom and return the first match in the list
    # remember that the lists are [match_idx, length/dist, reverse]
    # refer to docstring for decision making here
    if len(longest_protocore) == 1:
        match_idx = longest_protocore[0][0]

    elif len(central_protocore) > 0:
        match_idx = central_protocore[0][0]

    elif len(longest) == 1:
        match_idx = longest[0][0]

    elif len(central) > 0:
        match_idx = central[0][0]

    else:
        # this should never happen
        raise RuntimeError("No match found in LCS.")

    # match bounds on domain level will be inclusive to ensure a consistent exclusive
    # stop on cds level
    relevant_match = matching_block_dirs[match_idx]
    a_start = relevant_match[0]
    a_stop = relevant_match[0] + relevant_match[2] - 1  # inclusive stop
    b_start = relevant_match[1]
    b_stop = relevant_match[1] + relevant_match[2] - 1  # inclusive stop
    reverse = relevant_match[3]

    a_cds_start = a_domain_cds_idx[a_start]
    # cds stop may be end of cds
    if a_stop == len(a_domains) - 1:
        a_cds_stop = len(a_cds)
    else:
        a_cds_stop = a_domain_cds_idx[a_stop] + 1  # exclusive stop

    # fix b start and stop if in reverse. this means first flipping the domain indexes, getting the
    # cds index, and then flipping the cds index. yay.
    if reverse:
        old_start = b_start
        old_stop = b_stop
        b_start = len(b_domains) - b_stop - 1
        b_stop = len(b_domains) - old_start - 1  # inclusive stop

    b_cds_start = b_domain_cds_idx[b_start]

    # cds stop may be end of cds
    if b_stop == len(b_domains) - 1:
        b_cds_stop = len(b_cds)
    else:
        b_cds_stop = b_domain_cds_idx[b_stop] + 1  # exclusive stop

    if reverse:
        old_cds_start = b_cds_start
        b_cds_start = len(b_cds) - b_cds_stop
        b_cds_stop = len(b_cds) - old_cds_start  # exclusive stop

    # final check: it could happen that the start and stop of the domain LCS is in the
    # same CDS. in this case, the stop needs to be incremented by one
    if a_cds_start == a_cds_stop:
        a_cds_stop += 1

    if b_cds_start == b_cds_stop:
        b_cds_stop += 1

    # revert inclusive domain stop and index flipping before returning
    if reverse:
        b_start = old_start
        b_stop = old_stop
    a_stop += 1
    b_stop += 1

    return (
        a_start,
        a_stop,
        b_start,
        b_stop,
        a_cds_start,
        a_cds_stop,
        b_cds_start,
        b_cds_stop,
        reverse,
    )


def construct_missing_global_lcs(records: list[bs_genbank.BGCRecord], exemplar: int):
    """calculate missing lcs between family exemplar and family members

    Only relevant in GLOBAL alignment mode, as no lcs were calculated during comparisons

    Args:
        records (list[bs_genbank.BGCRecord]): family member records
        exemplar (int): idx of family exemplar record
    """
    record_db_dict = {record._db_id: record for record in records}

    # construct a list of member db ids, and separately exemplar db id
    db_ids = list(record_db_dict.keys())
    exemplar_id = db_ids.pop(exemplar)

    if DB.metadata is None:
        raise RuntimeError("DB metadata is None!")

    distance_table = DB.metadata.tables["distance"]

    # find exemplar->member pairs that do not have an lcs calculated yet
    missing_lcs_query = (
        select(distance_table.c.record_a_id, distance_table.c.record_b_id)
        .where(
            or_(
                and_(
                    distance_table.c.record_a_id == exemplar_id,
                    distance_table.c.record_b_id.in_(db_ids),
                ),
                and_(
                    distance_table.c.record_a_id.in_(db_ids),
                    distance_table.c.record_b_id == exemplar_id,
                ),
            )
        )
        .where(distance_table.c.lcs_domain_a_start == 0)
        .where(distance_table.c.lcs_domain_b_start == 0)
        .where(distance_table.c.reverse == 0)
    )

    missing_lcs_pairs = DB.execute(missing_lcs_query).fetchall()

    for missing_pair in missing_lcs_pairs:
        rec_a_id, rec_b_id = missing_pair
        pair = bs_comparison.RecordPair(
            record_db_dict[rec_a_id], record_db_dict[rec_b_id]
        )
        if isinstance(pair.record_a, bs_genbank.ProtoCluster) and isinstance(
            pair.record_b, bs_genbank.ProtoCluster
        ):
            lcs_data = find_domain_lcs_protocluster(pair)
        else:
            lcs_data = find_domain_lcs_region(pair)

        DB.execute(
            update(distance_table)
            .where(distance_table.c.record_a_id == rec_a_id)
            .where(distance_table.c.record_b_id == rec_b_id)
            .values(
                lcs_domain_a_start=lcs_data[0],
                lcs_domain_b_start=lcs_data[2],
                reverse=lcs_data[8],
            )
            .compile()
        )
