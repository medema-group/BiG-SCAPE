"""Module to handle parsing features from SeqIO records in GBK files"""

# from python


# from dependencies


# from other modules
import big_scape.enums as bs_enums
from big_scape.dereplicating.gbk_components import CDS
from big_scape.dereplicating.gbk_components import GBK


# from this module


antismash_parser_functions = {
    bs_enums.FEATURE_TYPE.CDS: CDS.parse
}
#     proto_core
#     protocluster
#     cand_cluster
#     region
#     cluster


minimal_parser_functions = {
     bs_enums.FEATURE_TYPE.CDS: CDS.parse
}


def get_parser_functions(run_mode: bs_enums.RUN_MODE) -> dict:
    """Function that chooses and returns correct parsing dict based on run mode and antiSMASH version

    Args:
        run_mode (bs_enums.RUN_MODE): cluster, query, or dereplicate
        as_version (str): antiSMASH version
        force_gbk (bool): force-gbk run param

    Returns:
        dict[str: function]: dict containing parsing functions for each seqIO feature type
    """

    if run_mode == bs_enums.RUN_MODE.DEREPLICATE:
        return minimal_parser_functions

    else:
        return antismash_parser_functions


def validate_cds_component(gbk: GBK) -> bool:
    """Checks whether the GBK contains at least 1 CDS component

    Returns:
        bool: True if the GBK contains at least 1 CDS component, False otherwise
    """

    return bs_enums.COMPONENTS.CDS in gbk.components
