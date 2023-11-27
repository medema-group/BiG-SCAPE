"""Module to generate newick GCF trees"""

# from python
import subprocess
import logging
import numpy as np
from typing import TextIO
from Bio import Phylo
from pathlib import Path
from collections import defaultdict
from scipy.optimize import linear_sum_assignment

# from other modules
from big_scape.genbank import BGCRecord


def generate_newick_tree(
    records: list[BGCRecord],
    exemplar: int,
    family_members: list[int],
    family_name: str,
    output_path: Path,
) -> str:
    """Generate newick formatted tree for each GCF

    Args:
        records (list[BGCRecord]): list of records within GCF
        exemplar (int): index of exemplar to use during alignment
        family_members (list[int]): list of bgc ids in one family
        family_name (str): name of current family
        output_path (Path): folder to store alignments and trees

    Returns:
        str: Correctly formatted newick tree
    """
    # no need for alignment
    if len(records) < 3:
        tree = f"({','.join([str(bgc_id)+':0.0' for bgc_id in family_members])}):0.01;"
    else:
        algn = generate_gcf_alignment(records, exemplar, family_members)
        algn_filename = output_path / Path(family_name + "_alignment.fasta")
        with open(algn_filename, "w") as out_algn:
            out_algn.write(algn)
        tree_filename = output_path / Path(family_name + ".newick")
        with open(tree_filename, "w") as out_newick:
            run_fasttree(algn_filename, out_newick)
        tree = process_newick_tree(tree_filename)
    return tree


def run_fasttree(algn_file: Path, out_file: TextIO):
    """Generate FastTree newick GCF tree

    Args:
        algn_file (Path): Path to alignment file
        out_file (TextIO): Opened output file object
    """

    result = subprocess.run(
        ["fasttree", "-quiet", "-nopr", algn_file],
        capture_output=True,
        shell=False,
    )
    stdout = result.stdout.decode("utf-8")
    out_file.write(stdout)

    if result.stderr:
        stderr = result.stderr.decode("utf-8")
        stderr = stderr.replace("\n", " ")
        logging.debug(f"FastTree says: {stderr}")


def process_newick_tree(tree_file: Path) -> str:
    """Process newick tree file format

    Args:
        tree_file (Path): Path to tree file

    Returns:
        str: processed newick tree
    """
    if not tree_file.exists():
        logging.error("Failed to create newick tree")
        raise FileNotFoundError()
    with open(tree_file, "r") as newick_file:
        try:
            tree = Phylo.read(newick_file, "newick")
        except ValueError as e:
            logging.warning("Error encountered while reading newick tree: ", str(e))
            return ""
        else:
            try:
                tree.root_at_midpoint()
            except UnboundLocalError:
                # Noticed this could happen if the sequences are exactly
                # the same and all distances == 0
                logging.debug("Unable to root at midpoint")
            return tree.format("newick")


def generate_gcf_alignment(
    records: list[BGCRecord], exemplar: int, family_members: list[int]
) -> str:
    """Generate protein domain alignment for records in GCF

    Args:
        records (list[BGCRecord]): Records within one GCF to align
        exemplar (int): Index of exemplar to use during alignment
        family_members (list[int]): list of bgc ids in one family

    Returns:
        str: alignment of GCF based on protein domain
        TODO: refactor
    """
    record_ids = list(range(len(records)))

    # collect present domains for each GCF member
    domain_sets = {}
    # count the frequency of occurrence of each domain (excluding copies)
    frequency_table: dict[str, int] = defaultdict(int)
    for idx, record in enumerate(records):
        domain_sets[idx] = set([domain.domain for domain in record.get_hsps()])
        for domain in domain_sets[idx]:
            frequency_table[domain] += 1

    # Find the set of [(tree domains)]. They should 1) be in the exemplar
    # and 2) appear with the most frequency. Iterate over the different
    # frequencies (descending) until set is not empty
    tree_domains: set[str] = set()
    frequencies = sorted(set(frequency_table.values()), reverse=True)

    # first try with domain(s) with max frequency, even if it's just one
    f = 0
    while len(tree_domains) == 0 and f < len(frequencies):
        for domain in frequency_table:
            if (
                frequency_table[domain] == frequencies[f]
                and domain in domain_sets[exemplar]
            ):
                tree_domains.add(domain)
        f += 1
    if len(tree_domains) == 1:
        logging.debug(
            "core shared domains for GCF {} consists of a single domain ({})".format(
                exemplar, [x for x in tree_domains][0]
            )
        )

    alignments: dict[int, str] = {}
    alignments[exemplar] = ""

    # store number of missed domains in bgc wrt exemplar
    missed_domains: dict[int, int] = {}
    record_ids.remove(exemplar)
    for record_idx in record_ids:
        alignments[record_idx] = ""
        missed_domains[record_idx] = 0

    match_dict: dict[int, int] = {}
    for domain in tree_domains:
        specific_domain_list_a = [
            hsp for hsp in records[exemplar].get_hsps() if hsp.domain == domain
        ]
        num_copies_a = len(specific_domain_list_a)
        for hsp in specific_domain_list_a:
            if hsp.alignment is not None:
                alignments[exemplar] += hsp.alignment.align_string
                seq_length = len(hsp.alignment.align_string)  # TODO: find better spot

        for bgc in alignments:
            match_dict.clear()
            if bgc == exemplar:
                pass
            elif domain not in domain_sets[bgc]:
                missed_domains[bgc] += 1
                alignments[bgc] += "-" * seq_length * num_copies_a
            else:
                specific_domain_list_b = [
                    hsp for hsp in records[bgc].get_hsps() if hsp.domain == domain
                ]
                num_copies_b = len(specific_domain_list_b)
                dist_matrix: np.ndarray = np.ndarray((num_copies_a, num_copies_b))

                for domsa in range(num_copies_a):
                    for domsb in range(num_copies_b):
                        hsp_a = specific_domain_list_a[domsa]
                        hsp_b = specific_domain_list_b[domsb]

                        if hsp_a.alignment is None or hsp_b.alignment is None:
                            logging.error(
                                "Trying to compare unaligned domains", hsp_a, hsp_b
                            )
                            raise AttributeError()
                        aligned_seq_a = hsp_a.alignment.align_string
                        aligned_seq_b = hsp_b.alignment.align_string

                        matches = 0
                        gaps = 0

                        for position in range(seq_length):
                            if aligned_seq_a[position] == aligned_seq_b[position]:
                                if aligned_seq_a[position] != "-":
                                    matches += 1
                                else:
                                    gaps += 1

                        dist_matrix[domsa][domsb] = 1 - (matches / (seq_length - gaps))

                best_indexes = linear_sum_assignment(dist_matrix)

                # at this point is not ensured that we have the same order
                # for the exemplar's copies (rows in BestIndexes)
                # ideally they should go from 0-numcopies. Better make sure

                for x in range(len(best_indexes[0])):
                    match_dict[best_indexes[0][x]] = best_indexes[1][x]

                for copy in range(num_copies_a):
                    try:
                        hsp_b = specific_domain_list_b[match_dict[copy]]
                        if hsp_b.alignment is None:
                            logging.error("Encountered unaligned domain", hsp_b)
                            raise AttributeError()
                        alignments[bgc] += hsp_b.alignment.align_string
                    except KeyError:
                        # This means that this copy of exemplar did not
                        # have a match in bgc (i.e. bgc has less copies
                        # of this domain than exemplar)
                        alignments[bgc] += "-" * seq_length

    # if a bgc is missing all tree domains, remove it from the tree
    delete_bgc: set[int] = set()
    for bgc in alignments:
        if bgc != exemplar and missed_domains[bgc] == len(tree_domains):
            delete_bgc.add(bgc)
    for bgc in delete_bgc:
        del alignments[bgc]

    algn_string = f">{family_members[exemplar]}\n{alignments[exemplar]}\n"
    for bgc in alignments:
        if bgc != exemplar:
            algn_string += f">{family_members[bgc]}\n{alignments[bgc]}\n"
    return algn_string
