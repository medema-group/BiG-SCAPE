"""Module to handle output generation for BiG-SCAPE dereplicate"""

# from python
import os
from pathlib import Path
import logging

# from dependencies

# from other modules
from big_scape.dereplicating.networking import Network
from big_scape.dereplicating.gbk_components.gbk import GBK


def write_output(network: Network, run: dict, gbk_list: list[GBK]) -> None:
    """Write output files for BiG-SCAPE dereplicate
    This function creates a directory for representative BGCs, creates symlinks
    for the representative BGCs, and writes a clustering tsv file that maps
    each BGC to its representative BGC.

    Args:
        network (Network): Network object containing connected components
        run (dict): run params
        gbk_list (list[GBK]): list of input GBK objects
    """

    output_dir: Path = run["output_dir"]
    input_dir: Path = run["input_dir"]

    rep_bgc_dir: Path = output_dir / "representative_clusters"

    # TODO: standardize existing output directory decisions
    # i.e. when to re-use and when to re-write.
    if not rep_bgc_dir.exists():
        logging.info("Creating representative BGCs directory: %s", rep_bgc_dir)
        rep_bgc_dir.mkdir(parents=False, exist_ok=True)
    else:
        logging.info("Representative BGCs directory already exists: %s", rep_bgc_dir)

    # create a mapping from the name of the GBK to the GBK object
    # so we can get the actual GBK objects with their full paths
    name_to_gbk = {}
    for gbk in gbk_list:
        gbk_name = gbk.name
        name_to_gbk[gbk_name] = gbk

    # create clustering tsv file
    clustering_tsv_path = output_dir / "clustering.tsv"
    logging.info("Writing clustering tsv file: %s", clustering_tsv_path)

    with clustering_tsv_path.open("w") as clustering_file:
        # write header
        clustering_file.write("gbk\trepresentative_gbk\n")

        for rep_name, cc in network.connected_components.items():
            rep_gbk: GBK = name_to_gbk.get(rep_name)

            # create a symlink for the representative BGC
            create_symlink(rep_gbk, rep_bgc_dir, input_dir)

            # write to the clustering tsv file
            for node_name in cc:
                node_gbk: GBK = name_to_gbk.get(node_name)
                clustering_file.write(
                    f"{node_gbk.path}\t{rep_gbk.path}\n"
                )


def create_symlink(rep_gbk: GBK, rep_bgc_dir: Path, input_dir: Path) -> None:
    """Create a symlink for the representative BGC in the output directory"""

    # absolute path to the representative GBK
    rep_gbk_path: Path = rep_gbk.path
    # path from input dir onwards
    rep_gbk_path_relative = rep_gbk_path.relative_to(input_dir)

    # create the symlink in the representative BGCs directory
    rep_gbk_symlink_path = rep_bgc_dir / rep_gbk_path_relative
    # ensure the parent directory exists
    rep_gbk_symlink_path.parent.mkdir(parents=True, exist_ok=True)

    if not rep_gbk_symlink_path.exists():
        os.symlink(rep_gbk_path, rep_gbk_symlink_path)
