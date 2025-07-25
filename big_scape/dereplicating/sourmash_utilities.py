"""Module to handle calculating distances between BGCs using sourmash and other required utilities"""

# from python
from pathlib import Path
import logging
import os
import subprocess
from subprocess import CalledProcessError

# from other modules
from big_scape.cli.config import BigscapeConfig
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

    sourmash_dir: Path = output_dir / "sourmash"

    # crash the run if sourmash_dir already exists
    if sourmash_dir.is_dir():
        logging.warning(
            "Sourmash output folder already exists: %s",
            sourmash_dir,
        )

    else:
        sourmash_dir.mkdir(parents=False, exist_ok=False)
        logging.info("Created sourmash output folder %s", sourmash_dir)

    CDS_fasta_dir: Path = Path(os.path.join(sourmash_dir, "CDS_fastas"))
    if not CDS_fasta_dir.is_dir():  # this should always be the case
        CDS_fasta_dir.mkdir(parents=False, exist_ok=False)
        logging.info("Created sourmash CDS fasta folder %s", CDS_fasta_dir)

    # create manysketch csv file
    manysketch_csv_path = Path(os.path.join(sourmash_dir, "manysketch.csv"))
    if manysketch_csv_path.is_file():
        logging.warning(
            "Manysketch csv file %s already exists, overwriting.",
            manysketch_csv_path,
        )
    with open(manysketch_csv_path, "w") as manysketch_csv_file:
        manysketch_csv_file.write("name,genome_filename,protein_filename\n")
        logging.info("Created manysketch csv file %s", manysketch_csv_path)

    logging.info("Concatenating CDSs and writing FASTA files to %s", CDS_fasta_dir)
    for gbk in gbk_list:
        CDS.concatenate_cds(gbk)
        fasta_file_path = get_fasta_path(gbk.name, CDS_fasta_dir)
        write_concat_cds_fasta(gbk, fasta_file_path)
        update_sourmash_fasta_manyketch_csv(
            manysketch_csv_path, gbk.name, fasta_file_path
        )

    return (sourmash_dir, CDS_fasta_dir, manysketch_csv_path)


def get_fasta_path(gbk_name: str, sourmash_dir: Path) -> tuple:
    """Get the names and paths of the fasta files

    Args:
        gbk_name (GBK): GBK object
        sourmash_dir (Path): path to the sourmash directory
    Returns:
        tuple: (gbk_name, fasta_file_path)
    """

    fasta_file_name = f"{gbk_name}.fasta"
    fasta_file_path = Path(os.path.join(sourmash_dir, fasta_file_name))

    return fasta_file_path


def update_sourmash_fasta_manyketch_csv(
    manysketch_csv_path: Path, gbk_name: str, fasta_file_path: Path
) -> None:
    """Update the manysketch csv file with the new fasta file

    Args:
        manysketch_csv_path (Path): _path to the manysketch csv file
        gbk_name (str): gbk name
        fasta_file_path (Path): path to the fasta file

    Returns:
        None
    """

    # write to manysketch csv file
    with open(manysketch_csv_path, "a") as manysketch_csv_file:
        # name, genome_filename, protein_filename
        # genome_filename is empty since we have only protein
        manysketch_csv_file.write(f"{gbk_name},,{fasta_file_path}\n")

    return None


def write_concat_cds_fasta(gbk: GBK, fasta_file_path: Path) -> None:
    """Write the concatenated CDS FASTA file for the given GBK object

    Args:
        gbk (GBK): GBK object
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
        fasta_file.write(f">{gbk.name}\n")
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
            "Sourmash sketch file %s already exists, overwriting.",
            sketch_file_path,
        )

    kmer_size = BigscapeConfig.KMER_SIZE
    scaled = BigscapeConfig.SCALED
    sketch_cmd = [
        "sourmash scripts manysketch",
        "-o",
        str(sketch_file_path),
        "-p",
        f"protein,k={kmer_size},scaled={scaled},noabund",
        "--singleton",  # one sketch per protein sequence
        "-c",
        str(run_dict["cores"]),
        str(manysketch_csv_path),
    ]

    sketch_cmd = " ".join(sketch_cmd)

    sketch_log_path = Path(os.path.join(sourmash_dir, "sketch.log"))
    with open(sketch_log_path, "w") as sketch_log_file:
        logging.info("Running sourmash sketch")
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
                "Sourmash sketch failed with error code %d: %s",
                e.returncode,
                e.stderr,
            )
            sketch_log_file.write(e.stderr)
            raise e
    logging.info("Sourmash sketch completed successfully")

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
            "Sourmash pairwise distance file %s already exists, overwriting.",
            pairwise_file_path,
        )

    kmer_size = BigscapeConfig.KMER_SIZE
    scaled = BigscapeConfig.SCALED
    pairwise_cmd = [
        "sourmash scripts pairwise",
        "-o",
        str(pairwise_file_path),
        "-k",  # kmer size
        str(kmer_size),
        "-s",  # scaling factor
        str(scaled),
        "-m",
        "protein",  # sequence type
        "--write-all",  # write self comparisons for all sketches
        "-A",  # ignore containment threshold and output all comparisons
        "-c",
        str(run_dict["cores"]),
        str(sketch_file_path),
    ]

    pairwise_cmd = " ".join(pairwise_cmd)

    pairwise_log_path = Path(os.path.join(sourmash_dir, "pairwise.log"))
    with open(pairwise_log_path, "w") as pairwise_log_file:
        logging.info("Running sourmash pairwise")
        logging.debug("Sourmash pairwise command: %s", pairwise_cmd)
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
                "Sourmash pairwise failed with error code %d: %s",
                e.returncode,
                e.stderr,
            )
            pairwise_log_file.write(e.stderr)
            raise e
    logging.info("Sourmash pairwise completed successfully")

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

    edges = set()
    nodes = set()

    with open(pairwise_file_path, "r") as pairwise_file:
        for line in pairwise_file:
            if line.startswith("query"):
                continue

            parts = line.strip().split(",")
            nodeA = parts[0]
            nodeB = parts[2]
            distance = parts[6]  # jaccard similarity

            # skip self-comparisons
            # we dont want to rely on this being present in the file
            if nodeA == nodeB:
                continue

            # add nodes to the set of nodes
            nodes.update((nodeA, nodeB))

            edge = Edge(nodeA, nodeB, float(distance))
            if edge.jaccard_similarity >= cutoff:
                edges.add(edge)

    logging.info("Parsed %d edges from sourmash results (similarity threshold: %.2f)", len(edges), cutoff)

    edges = list(edges)
    nodes = list(nodes)

    return edges, nodes
