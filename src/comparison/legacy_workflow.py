"""Contains methods to run the legacy comparison workflow on a bin of BGC pairs"""

# from python
import logging

# from other modules
from src.distances import calc_jaccard_pair, calc_ai_pair, calc_dss_pair_legacy
from src.network import BSNetwork

# from this module
from .legacy_extend import (
    legacy_needs_extend,
    expand_glocal,
    check_expand,
    reset_expansion,
)
from .binning import BGCBin


def create_bin_network_edges_old(bin: BGCBin, network: BSNetwork, alignment_mode: str):
    for pair in bin.pairs(legacy_sorting=True):
        # calculate jaccard for the full sets. if this is 0, there are no shared domains
        # important not to cache here otherwise we are using the full range again later
        jaccard = calc_jaccard_pair(pair, cache=False)

        if jaccard == 0.0:
            network.add_edge(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

        logging.debug("JC: %f", jaccard)

        # we record the LCS starts and stops here in case we need to reset them later
        (
            lcs_a_start,
            lcs_a_stop,
            lcs_b_start,
            lcs_b_stop,
        ) = pair.comparable_region.find_lcs()
        pair.comparable_region.log_comparable_region("LCS")

        if legacy_needs_extend(pair, alignment_mode):
            expand_glocal(pair.comparable_region)

            if check_expand(pair.comparable_region):
                pair.comparable_region.log_comparable_region("GLOCAL")

                jaccard = calc_jaccard_pair(pair)

                if jaccard == 0.0:
                    network.add_edge(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
                    continue
            else:
                reset_expansion(pair.comparable_region)
        else:
            reset_expansion(pair.comparable_region)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

        # mix
        distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

        logging.debug(
            "JC: %f, AI: %f, DSS: %f, SCORE: %f", jaccard, adjacency, dss, distance
        )

        network.add_edge(pair, jc=jaccard, ai=adjacency, dss=dss, dist=distance)


def create_bin_network_edges(bin: BGCBin, network: BSNetwork, alignment_mode: str):
    logging.info("Calculating Jaccard for %d pairs", bin.num_pairs())
    for pair in bin.pairs(legacy_sorting=True):
        # calculate jaccard for the full sets. if this is 0, there are no shared domains
        # important not to cache here otherwise we are using the full range again later
        jaccard = calc_jaccard_pair(pair, cache=False)

        if jaccard == 0.0:
            network.add_edge(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

    logging.info(
        "Performing LCS for %d pairs", bin.num_pairs() - network.graph.number_of_edges()
    )
    pairs_need_expand = []
    pairs_no_expand = []
    for pair in bin.pairs(legacy_sorting=True):
        if pair in network:
            continue

        logging.debug("JC: %f", jaccard)

        pair.comparable_region.find_lcs()
        pair.comparable_region.log_comparable_region("LCS")

        if legacy_needs_extend(pair, alignment_mode):
            pairs_need_expand.append(pair)
            continue

        reset_expansion(pair.comparable_region)
        pairs_no_expand.append(pair)

    logging.info("Expanding regions for %d pairs", len(pairs_need_expand))
    expanded_pairs = []
    for pair in pairs_need_expand:
        expand_glocal(pair.comparable_region)

        if not check_expand(pair.comparable_region):
            reset_expansion(pair.comparable_region)
            pairs_no_expand.append(pair)
            continue

        pair.comparable_region.log_comparable_region("GLOCAL")

        jaccard = calc_jaccard_pair(pair)

        if jaccard == 0.0:
            network.add_edge(pair, jc=0.0, ai=0.0, dss=0.0, dist=1.0)
            continue

        expanded_pairs.append(pair)

    logging.info(
        "Calculating score for %d pairs that were reset or did not need expansion",
        len(pairs_no_expand),
    )
    for pair in pairs_no_expand:
        jaccard = calc_jaccard_pair(pair)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

        # mix
        distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

        logging.debug(
            "JC: %f, AI: %f, DSS: %f, SCORE: %f", jaccard, adjacency, dss, distance
        )

        network.add_edge(pair, jc=jaccard, ai=adjacency, dss=dss, dist=distance)

    logging.info(
        "Calculating score for %d pairs that were expanded",
        len(expanded_pairs),
    )
    for pair in expanded_pairs:
        jaccard = calc_jaccard_pair(pair)

        adjacency = calc_ai_pair(pair)
        # mix anchor boost = 2.0
        dss = calc_dss_pair_legacy(pair, anchor_boost=2.0)

        # mix
        distance = 1 - (0.2 * jaccard) - (0.05 * adjacency) - (0.75 * dss)

        logging.debug(
            "JC: %f, AI: %f, DSS: %f, SCORE: %f", jaccard, adjacency, dss, distance
        )

        network.add_edge(pair, jc=jaccard, ai=adjacency, dss=dss, dist=distance)
