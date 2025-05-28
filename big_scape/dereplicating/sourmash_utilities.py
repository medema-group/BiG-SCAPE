"""Module to handle calculating distances between BGCs using sourmash and other required utilities"""

# from python
from pathlib import Path
import logging
import os
import subprocess
from subprocess import CalledProcessError

# from other modules
from big_scape.dereplicating.gbk_components.gbk import GBK
from big_scape.dereplicating.gbk_components.cds import CDS
import big_scape.enums as bs_enums
from big_scape.dereplicating.networking import Edge


def make_sourmash_input(gbk_list: list[GBK], run: dict) -> Path:
    """Create the sourmash input files
    Args:
        gbk_list (list[GBK]): list of GBK objects
        run (dict): run dictionary containing the input and output directories
    Returns:
        sourmash_fasta_folder (Path): path to the sourmash fasta folder
    """

    output_dir: Path = run["output_dir"]

    if not output_dir.is_dir():
        output_dir.mkdir(parents=False, exist_ok=False)

    sourmash_dir: Path = output_dir/"sourmash"

    # crash the run if sourmash_dir already exists
    if sourmash_dir.is_dir():
        logging.warning(
            "Sourmash output folder %s already exists, exiting.",
            sourmash_dir,
        )
        # TODO: consider desired behaviour if sourmash_dir already exists
        # re-write files? crash the run? something in between?
        # raise FileExistsError(
        #     f"Sourmash output folder {sourmash_dir} already exists, exiting."
        # )

    else:
        sourmash_dir.mkdir(parents=False, exist_ok=False)
        logging.info("Created sourmash output folder %s", sourmash_dir)

    CDS_fasta_dir: Path = Path(os.path.join(sourmash_dir, "CDS_fastas"))
    if not CDS_fasta_dir.is_dir():  # this should always be the case
        CDS_fasta_dir.mkdir(parents=False, exist_ok=False)
        logging.info("Created sourmash CDS fasta output folder %s", CDS_fasta_dir)

    # create manysketch csv file
    manysketch_csv_path = Path(os.path.join(sourmash_dir, "manysketch.csv"))
    with open(manysketch_csv_path, "w") as manysketch_csv_file:
        manysketch_csv_file.write("name,genome_filename,protein_filename\n")
        logging.info("Created manysketch csv file %s", manysketch_csv_path)

    logging.info("Concatenating CDSs and writing FASTA files to %s", CDS_fasta_dir)
    for gbk in gbk_list:
        CDS.concatenate_cds(gbk)
        gbk_name, fasta_file_path = get_fasta_names_paths(gbk, run, CDS_fasta_dir)
        write_concat_cds_fasta(gbk, gbk_name, fasta_file_path)
        update_sourmash_fasta_manyketch_csv(
            manysketch_csv_path, gbk_name, gbk.hash, fasta_file_path
        )

    return (sourmash_dir, CDS_fasta_dir, manysketch_csv_path)


def get_fasta_names_paths(gbk: GBK, run: dict, sourmash_dir: Path) -> tuple:
    """Get the names and paths of the fasta files

    Args:
        gbk (GBK): GBK object
        run (dict): run dictionary containing the input and output directories
        sourmash_dir (Path): path to the sourmash directory
    Returns:
        tuple: (gbk_name, fasta_file_path)
    """

    input_dir = run["input_dir"]

    # get path from input dir onwards
    gbk_path = gbk.path.relative_to(input_dir)

    # TODO: consider if this is the best way to name the fasta files
    # gbk name will be: path from input dir onwards
    # and the hash of the gbk
    gbk_name = ".".join(gbk_path.parts)

    fasta_file_name = f"{gbk_name}.{gbk.hash}.fasta"
    fasta_file_path = Path(os.path.join(sourmash_dir, fasta_file_name))

    return (gbk_name, fasta_file_path)


def update_sourmash_fasta_manyketch_csv(
    manysketch_csv_path: Path, gbk_name: str, gbk_hash: str, fasta_file_path: Path
) -> None:
    """Update the manysketch csv file with the new fasta file

    Args:
        manysketch_csv_path (Path): _path to the manysketch csv file
        gbk_name (str): gbk name
        gbk_hash (str): gbk hash
        fasta_file_path (Path): path to the fasta file

    Returns:
        None
    """

    # write to manysketch csv file
    with open(manysketch_csv_path, "a") as manysketch_csv_file:
        # name, genome_filename, protein_filename
        # genome_filename is empty since we have only protein
        manysketch_csv_file.write(f"{gbk_name}.{gbk_hash},,{fasta_file_path}\n")

    return None


def write_concat_cds_fasta(gbk: GBK, gbk_name: str, fasta_file_path: Path) -> None:
    """Write the concatenated CDS FASTA file for the given GBK object

    Args:
        gbk (GBK): GBK object
        gbk_name (str): name of the GBK object
        fasta_file_path (Path): path to the FASTA file
    Returns:
        None
    """

    concat_cds = gbk.components[bs_enums.COMPONENTS.CONCAT_CDS]

    if concat_cds is None:
        logging.warning("No concat CDS for %s, skipping.", gbk.path)
        return None

    if fasta_file_path.is_file():  # this should never happen
        logging.warning(
            "Trying to write CDS sequence of %s to %s, but file already exists, skipping.",
            gbk.path,
            fasta_file_path,
        )
        return None

    seq = "\n".join(
        str(concat_cds)[i: i + 80] for i in range(0, len(str(concat_cds)), 80)
    )

    with open(fasta_file_path, "w") as fasta_file:
        logging.debug(
            "Writing concatenated CDS sequence of %s to %s", gbk.path, fasta_file_path
        )
        fasta_file.write(f">{gbk_name}.{gbk.hash}\n")
        fasta_file.write(f"{seq}\n")

    return None


def run_sourmash_branchwater(
    run_dict: dict, sourmash_dir: Path, CDS_fastas_dir: Path, manysketch_csv_path: Path
) -> Path:
    """Run sourmash branchwater plugin

    Args:
        run_dict (dict): run dictionary
        sourmash_dir (Path): path to the sourmash directory
        CDS_fastas_dir (Path): path to the CDS fasta directory
        manysketch_csv_path (Path): path to the manysketch csv file
    Returns:
        sourmash_distances_file_path (Path): path to the sourmash distances file
    """

    sourmash_sketch_file_path = sourmash_sketch(
        run_dict, sourmash_dir, CDS_fastas_dir, manysketch_csv_path
    )

    sourmash_distances_file_path = sourmash_compare(
        run_dict, sourmash_dir, sourmash_sketch_file_path
    )

    return sourmash_distances_file_path


def sourmash_sketch(
    run_dict: dict, sourmash_dir: Path, CDS_fastas_dir: Path, manysketch_csv_path: Path
) -> Path:
    """Run sourmash sketch

    Args:
        run_dict (dict): run dictionary
        sourmash_dir (Path): path to the sourmash directory
        CDS_fastas_dir (Path): path to the CDS fasta directory
        manysketch_csv_path (Path): path to the manysketch csv file
    Returns:
        sketch_file_path (Path): path to the sourmash sketch file
    """

    sketch_file_path = Path(os.path.join(sourmash_dir, "sketch_proteins.zip"))

    if sketch_file_path.is_file():
        logging.warning(
            "Sourmash sketch file %s already exists, skipping.",
            sketch_file_path,
        )
        return sketch_file_path

    # TODO: consider adding sourmash params from config file
    # needs some discussion however on wether this is desireable
    sketch_cmd = [
        "sourmashscripts manysketch",
        "-o",
        str(sketch_file_path),
        "-p",
        "protein,k=10,scaled=200,noabund",  # defailt params for protein sketch
        "--singleton",
        "-c",
        str(run_dict["cores"]),
        str(manysketch_csv_path),
    ]

    sketch_cmd = " ".join(sketch_cmd)

    sketch_log_path = Path(os.path.join(sourmash_dir, "sketch.log"))
    with open(sketch_log_path, "w") as sketch_log_file:
        logging.info("Running branchwater sourmash sketch")
        logging.debug("Sourmash branchwater sketch command: %s", sketch_cmd)
        try:
            subprocess.run(
                sketch_cmd,
                shell=True,
                check=True,
                text=True,
                stdout=sketch_log_file,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
        except CalledProcessError as e:
            logging.error(
                "Sourmash branchwater sketch failed with error code %d: %s",
                e.returncode,
                e.stderr,
            )
            sketch_log_file.write(e.stderr)
            raise e
    logging.info("Sourmash branchwater sketch completed successfully")

    return sketch_file_path


def sourmash_compare(
    run_dict: dict, sourmash_dir: Path, sketch_file_path: Path
) -> None:
    """Run sourmash pairwise comparison

    Args:
        run_dict (dict): run dictionary
        sourmash_dir (Path): path to the sourmash directory
        sketch_file_path (Path): path to the sourmash sketch file
    Returns:
        pairwise_file_path (Path): path to the sourmash pairwise file
    """

    pairwise_file_path = Path(os.path.join(sourmash_dir, "pairwise_output.csv"))

    if pairwise_file_path.is_file():
        logging.warning(
            "Sourmash pairwise distance file %s already exists, skipping.",
            pairwise_file_path,
        )
        return pairwise_file_path

    pairwise_cmd = [
        "sourmash scripts pairwise",
        "-o",
        str(pairwise_file_path),
        "-k",
        "10",
        "-s",
        "200",
        "-m",
        "protein",
        "--write-all",  # write self comparisons for all sketches
        "-A",  # ignore containment threshold and output all comparisons
        "-c",
        str(run_dict["cores"]),
        str(sketch_file_path),
    ]

    pairwise_cmd = " ".join(pairwise_cmd)

    pairwise_log_path = Path(os.path.join(sourmash_dir, "pairwise.log"))
    with open(pairwise_log_path, "w") as pairwise_log_file:
        logging.info("Running sourmash branchwater pairwise")
        logging.debug("Sourmash branchwater pairwise command: %s", pairwise_cmd)
        try:
            subprocess.run(
                pairwise_cmd,
                shell=True,
                check=True,
                text=True,
                stdout=pairwise_log_file,
                stderr=subprocess.PIPE,
                encoding="utf-8",
            )
        except CalledProcessError as e:
            logging.error(
                "Sourmash branchwater pairwise failed with error code %d: %s",
                e.returncode,
                e.stderr,
            )
            pairwise_log_file.write(e.stderr)
            raise e
    logging.info("Sourmash branchwater pairwise completed successfully")

    return pairwise_file_path


def parse_sourmash_results(pairwise_file_path: Path, cutoff: float) -> set[Edge]:
    """Parse the sourmash pairwise results, and return a list of edges

    Args:
        pairwise_file_path (Path): path to the sourmash pairwise file
        cutoff (float): jaccard similarity cutoff for edges
    Returns:
        edges (set[Edge]): list of edges
    """

    if not pairwise_file_path.is_file():
        logging.error("Sourmash pairwise file %s does not exist", pairwise_file_path)
        raise FileNotFoundError()

    edges = []

    with open(pairwise_file_path, "r") as pairwise_file:
        for line in pairwise_file:
            if line.startswith("query"):
                continue

            nodeA, _, nodeB, _, _, _, distance, _, _, _, _ = line.strip().split(",")
            if nodeA == nodeB:
                # we dont want to rely on this being present in the file
                continue
            edge = Edge(nodeA, nodeB, float(distance))
            if edge.jaccard_similarity >= cutoff:
                edges.append(edge)
            else:
                # TODO: consider if this is the best way to make sure
                # that all nodes are present in the network, i.e.
                # if the edge is below the cutoff, we still want to add
                # the nodes to the network
                singletonA = Edge(nodeA, nodeA, 1.0)
                edges.append(singletonA)

                singletonB = Edge(nodeB, nodeB, 1.0)
                edges.append(singletonB)

    logging.info("Parsed %d edges from sourmash results", len(edges))

    edges = set(edges)  # remove duplicates
    return edges

