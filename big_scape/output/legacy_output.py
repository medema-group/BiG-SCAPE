"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
from itertools import repeat
import json
import logging
from multiprocessing import Pool
from distutils import dir_util
from pathlib import Path
import click
from sqlalchemy import select, alias
from typing import Optional

# from other modules
from big_scape.data import DB
from big_scape.comparison import RecordPairGenerator
from big_scape.genbank import GBK, BGCRecord
from big_scape.trees import generate_newick_tree, save_trees
from big_scape.comparison import lcs, get_record_category

import big_scape.paths as bs_paths
import big_scape.enums as bs_enums


def copy_base_output_templates(output_dir: Path):
    """Copy the base output html/javascript templates to an output
    directory

    Args:
        output_dir (Path): main output directory
    """

    template_root = bs_paths.TEMPLATES_OUTPUT_DIR / Path("html_template")
    template_dir = template_root / "output"

    # copy html content
    dir_util.copy_tree(str(template_dir), str(output_dir), verbose=False)


def generate_pfams_js(output_dir: Path, pfam_info: list[tuple[str, str, str]]) -> None:
    """Generate the pfam.js file needed to show the correct PFAM domain colors

    Args:
        output_dir (Path): main output directory
        pfam_info (list[tuple[str]]): A list of tuples containing pfam information. Each tuple
        contains an accession, a name and a description
    """

    # gather color information
    pfam_colors_file_path = bs_paths.TEMPLATES_OUTPUT_DIR / Path("domain_colors.tsv")

    # make accession to color dictionary
    pfam_colors_dict = {}
    with open(pfam_colors_file_path, mode="r", encoding="utf8") as pfam_colors_file:
        for line in pfam_colors_file:
            line_parts = line.rstrip().split("\t")
            pfam_colors_dict[line_parts[0]] = line_parts[1]

    pfams_data = {}

    for info in pfam_info:
        accession, name, description = info

        # trim version number
        accession = accession[:7]

        # default to white
        color = pfam_colors_dict.get(accession, "255,255,255")

        pfams_data[accession] = {
            "col": color,
            "name": name,
            "desc": description,
        }

    pfams_path = output_dir / Path("html_content/js/pfams.js")

    with open(pfams_path, "w") as run_data_js:
        run_data_js.write(
            "var pfams={};\n".format(
                json.dumps(pfams_data, indent=4, separators=(",", ":"), sort_keys=True)
            )
        )


def generate_bs_families_members(
    cutoff: float, pair_generator: RecordPairGenerator, run_id: int
) -> dict[int, list[int]]:
    """Generate a dictionary where keys are family indexes and values are list of region
    indexes that belong to that family

    Args:
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator
        run_id: id of the current run

    Returns:
        dict[int, list[int]]: family to member regions index
    """
    # get a dictionary of node id to family id
    node_family = {}

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    # get all families from the database
    family_table = DB.metadata.tables["family"]
    bgc_families_table = DB.metadata.tables["bgc_record_family"]

    select_statement = (
        select(
            bgc_families_table.c.record_id,
            family_table.c.center_id,
            family_table.c.cutoff,
            family_table.c.bin_label,
        )
        .join(family_table, bgc_families_table.c.family_id == family_table.c.id)
        .where(bgc_families_table.c.record_id.in_(pair_generator.record_ids))
        .where(family_table.c.cutoff == cutoff)
        .where(family_table.c.bin_label == pair_generator.label)
        .where(family_table.c.run_id == run_id)
    )

    result = DB.execute(select_statement).fetchall()

    for row in result:
        node_family[row.record_id] = row.center_id

    families_members: dict[int, list[int]] = {}

    for idx, record in enumerate(pair_generator.source_records):
        record_id = record._db_id

        if record_id is None:
            raise AttributeError("Record id is None!")

        if record_id not in node_family:
            families_members[record_id] = [idx]
            continue

        family_id = node_family[record_id]

        if family_id not in families_members:
            families_members[family_id] = []

        families_members[family_id].append(idx)

    return families_members


def legacy_prepare_output(
    output_dir: Path, pfam_info: list[tuple[str, str, str]]
) -> None:
    """Prepare the base output files for the run. These are not specific to cutoffs or
    to pair_generators

    Args:
        output_dir (Path): output directory
        pfam_info (list[tuple(str, str, str)]): A list of tuples containing information
        about the PFAM entries used in this analysis. tuples correspond to the accession
        name and description of a pfam entry
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        logging.info("Skipping all (non-SQLite db related) output generation")
        return

    copy_base_output_templates(output_dir)

    generate_pfams_js(output_dir, pfam_info)

    # generate root output_files folder
    output_files_root = output_dir / "output_files"
    if not output_files_root.exists():
        output_files_root.mkdir(exist_ok=True)


def legacy_prepare_cutoff_output(run: dict, cutoff: float, gbks: list[GBK]) -> None:
    """Prepare output data for a given cutoff value

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        gbks (list[GBK]): list of gbks used in the analysis
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        return

    output_dir: Path = run["output_dir"]
    label = run["label"]

    # networks subfolders
    output_files_root = output_dir / "output_files"
    cutoff_files_path = output_files_root / f"{label}_c{cutoff}"
    cutoff_files_path.mkdir(exist_ok=True)


def legacy_prepare_bin_output(
    run: dict, cutoff: float, pair_generator: RecordPairGenerator
) -> None:
    """Prepare output data for a given pair_generator at a given cutoff value

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        return

    output_dir: Path = run["output_dir"]
    label = run["label"]

    # networks subfolders
    output_files_root = output_dir / "output_files"
    cutoff_files_path = output_files_root / f"{label}_c{cutoff}"
    pair_generator_files_path = cutoff_files_path / pair_generator.label
    pair_generator_files_path.mkdir(exist_ok=True)

    if click_context and click_context.obj["no_trees"]:
        return

    tree_path = pair_generator_files_path / Path("GCF_trees")
    tree_path.mkdir(exist_ok=True)


def legacy_generate_bin_output(
    run: dict,
    cutoff: float,
    pair_generator: RecordPairGenerator,
) -> None:
    """Generate the network data from a pair_generator from cutoff filtering and affinity
    propagation

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator
        network (BSNetwork): the network object for the pair_generator
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        return

    families_members = generate_bs_families_members(
        cutoff, pair_generator, run["run_id"]
    )
    write_clustering_file(run, cutoff, pair_generator)
    write_cutoff_network_file(run, cutoff, pair_generator)

    if click_context and click_context.obj["no_trees"]:
        return

    generate_newick_trees(run, cutoff, pair_generator, families_members)


def generate_newick_trees(
    run: dict,
    cutoff: float,
    pair_generator: RecordPairGenerator,
    families_members: dict[int, list[int]],
) -> None:
    """Generate newick trees for all families

    Args:
        run (dict): current run parameters
        cutoff (float): current bin cutoff
        pair_generator (RecordPairGenerator): current bin BGC pair generator
        families_members (dict[int, list[int]]): map of family center db id to family
            member indices in the pair generators source records
    """
    tree_path = (
        run["output_dir"]
        / "output_files"
        / f"{run['label']}_c{cutoff}"
        / pair_generator.label
        / "GCF_trees"
    )

    exemplars_db = []
    families_records = []
    exemplars_idx = []
    for exemplar, members in families_members.items():
        family_recs = [pair_generator.source_records[member] for member in members]
        for idx, record in enumerate(family_recs):
            if record._db_id == exemplar:
                exemplar_idx = idx
                break
        exemplars_db.append(exemplar)
        families_records.append(family_recs)
        exemplars_idx.append(exemplar_idx)

    if run["alignment_mode"] == bs_enums.ALIGNMENT_MODE.GLOBAL:
        for records, exemplar in zip(families_records, exemplars_idx):
            lcs.construct_missing_global_lcs(records, exemplar)

    pool = Pool(processes=run["cores"])
    trees = pool.starmap(
        generate_newick_tree,
        zip(families_records, exemplars_idx, exemplars_db, repeat(tree_path)),
    )
    save_trees(trees, exemplars_db, cutoff, pair_generator.label, run["run_id"])
    DB.commit()


def write_record_annotations_file(run, cutoff, all_bgc_records) -> None:
    """Writes the record annotations file to the output directory
    Formerly Network_Annotations_Full.tsv

    Args:
        run (dict): contains the run parameters
        cutoff (float): distance metric cutoff
        all_bgc_records (list): contains all bgc records in this run

    Raises:
        RuntimeError: if database not available
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        return

    run_id = run["run_id"]
    output_dir = run["output_dir"]
    label = run["label"]
    output_files_root = output_dir / "output_files"
    cutoff_path = output_files_root / f"{label}_c{cutoff}"
    record_annotations_path = cutoff_path / "record_annotations.tsv"

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")
    bgc_record_table = DB.metadata.tables["bgc_record"]
    gbk_table = DB.metadata.tables["gbk"]
    cc_table = DB.metadata.tables["connected_component"]

    record_categories = {}
    for record in all_bgc_records:
        record_categories[record._db_id] = get_record_category(record)

    select_statement = (
        select(
            gbk_table.c.path,
            gbk_table.c.organism,
            gbk_table.c.taxonomy,
            gbk_table.c.description,
            bgc_record_table.c.id,
            bgc_record_table.c.record_number,
            bgc_record_table.c.record_type,
            bgc_record_table.c.product,
        )
        .join(bgc_record_table, bgc_record_table.c.gbk_id == gbk_table.c.id)
        .join(cc_table, cc_table.c.record_id == bgc_record_table.c.id)
        .where(cc_table.c.cutoff == cutoff)
        .where(cc_table.c.run_id == run_id)
        .distinct()
    )

    record_data = DB.execute(select_statement).fetchall()

    with open(record_annotations_path, "w") as record_annotations_file:
        header = "\t".join(
            [
                "GBK",
                "Record_Type",
                "Record_Number",
                "Class",
                "Category",
                "Organism",
                "Taxonomy",
                "Description",
            ]
        )
        record_annotations_file.write(header + "\n")

        for record in record_data:
            (
                gbk_path,
                organism,
                taxonomy,
                description,
                rec_id,
                record_number,
                record_type,
                product,
            ) = record

            row = "\t".join(
                [
                    Path(gbk_path).stem,
                    record_type,
                    str(record_number),
                    product,
                    record_categories[rec_id],
                    organism,
                    taxonomy,
                    description,
                ]
            )
            record_annotations_file.write(row + "\n")

    return None


def write_clustering_file(run, cutoff, pair_generator) -> None:
    """Writes the clustering file to the output directory
    bin_clustering_cutoff.tsv

    Args:
        run (dict): run parameters
        cutoff (float): distance metric cutoff
        pair_generator (RecordPairGenerator): BGC record bin

    Raises:
        RuntimeError: database not found
    """

    output_dir = run["output_dir"]
    label = run["label"]
    bin_label = pair_generator.label
    run_id = run["run_id"]

    output_files_root = output_dir / "output_files"
    cutoff_path = output_files_root / f"{label}_c{cutoff}"
    pair_generator_path = cutoff_path / pair_generator.label
    clustering_file_path = pair_generator_path / f"{bin_label}_clustering_c{cutoff}.tsv"

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    gbk_table = DB.metadata.tables["gbk"]
    bgc_record_table = DB.metadata.tables["bgc_record"]
    family_table = DB.metadata.tables["family"]
    rec_fam_table = DB.metadata.tables["bgc_record_family"]
    cc_table = DB.metadata.tables["connected_component"]

    record_ids = pair_generator.record_ids
    select_statement = (
        select(
            gbk_table.c.path,
            bgc_record_table.c.record_type,
            bgc_record_table.c.record_number,
            cc_table.c.id,
            rec_fam_table.c.family_id,
        )
        .join(bgc_record_table, bgc_record_table.c.gbk_id == gbk_table.c.id)
        .join(rec_fam_table, bgc_record_table.c.id == rec_fam_table.c.record_id)
        .join(family_table, rec_fam_table.c.family_id == family_table.c.id)
        .join(cc_table, cc_table.c.record_id == bgc_record_table.c.id)
        .where(rec_fam_table.c.record_id.in_(record_ids))
        .where(family_table.c.cutoff == cutoff)
        .where(family_table.c.bin_label == bin_label)
        .where(family_table.c.run_id == run_id)
        .where(cc_table.c.cutoff == cutoff)
        .where(cc_table.c.bin_label == bin_label)
        .where(cc_table.c.run_id == run_id)
        .order_by(rec_fam_table.c.family_id)
    )

    record_data = DB.execute(select_statement).fetchall()

    with open(clustering_file_path, "w") as clustering_file:
        header = "\t".join(["GBK", "Record_Type", "Record_Number", "CC", "Family"])
        clustering_file.write(header + "\n")

        for record in record_data:
            gbk_path, record_type, record_number, cc_number, family = record

            row = "\t".join(
                [
                    Path(gbk_path).stem,
                    record_type,
                    str(record_number),
                    str(cc_number),
                    f"FAM_{family:0>5}",
                ]
            )
            clustering_file.write(row + "\n")

    return None


def write_cutoff_network_file(
    run: dict, cutoff: float, pair_generator: RecordPairGenerator
) -> None:
    """Writes the cutoff network file to the output directory
    i.e. edge list for a given bin with edges above the cutoff

    Args:
        run (dict): run parameters
        cutoff (float): distance cutoff
        pair_generator (RecordPairGenerator): bin with BGC records

    """

    output_dir = run["output_dir"]
    label = run["label"]
    bin_label = pair_generator.label

    output_files_root = output_dir / "output_files"
    cutoff_path = output_files_root / f"{label}_c{cutoff}"
    pair_generator_path = cutoff_path / pair_generator.label
    cutoff_network_file_path = pair_generator_path / f"{bin_label}_c{cutoff}.network"

    write_network_file(
        run,
        pair_generator.source_records,
        cutoff_network_file_path,
        pair_generator.weights,
        cutoff,
    )


def write_full_network_file(run: dict, all_bgc_records: list[BGCRecord]) -> None:
    """Writes the full network file to the output directory,
    i.e. all edges from db that have both records in the run,
    and for all weights relevant to the run

    Args:
        run (dict): run parameters
        all_bgc_records (list): BGC records in this run
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["db_only_output"]:
        return

    output_dir = run["output_dir"]
    label = run["label"]

    output_files_root = output_dir / "output_files"
    full_network_file_path = output_files_root / f"{label}_full.network"

    write_network_file(run, all_bgc_records, full_network_file_path)


def write_network_file(
    run: dict,
    bgc_records: list[BGCRecord],
    output_path: Path,
    weight: Optional[str] = None,
    cutoff: Optional[float] = None,
) -> None:
    """Write all edges to a file for all pairs in a list of records

    If no weight is given, returns edges for all relevant weights in this run.
    If no cutoff is given, returns all relevant edges regardless of distance

    Args:
        run (dict): run parameters
        bgc_records (list[BGCRecord]): list of records to collect edges for
        output_path (Path): output network file path
        weight (str, optional): bin weight to collect edges for. Defaults to None.
        cutoff (float, optional): distance cutoff for returned edges. Defaults to None.
    """
    if weight is None:
        legacy_weights = [
            "PKSI",
            "PKSother",
            "NRPS",
            "RiPPs",
            "saccharides",
            "terpene",
            "PKS-NRP_Hybrids",
            "other",
        ]
        incl_weights = ["mix"]

        if not run["mix"]:
            incl_weights.remove("mix")
        if run["legacy_weights"]:
            incl_weights.extend(legacy_weights)
    else:
        incl_weights = [weight]

    aln_mode = run["alignment_mode"].name
    ext_strat = run["extend_strategy"].name

    record_ids = [record._db_id for record in bgc_records]

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    distance_table = DB.metadata.tables["distance"]
    gbk_table = DB.metadata.tables["gbk"]
    bgc_record_table = DB.metadata.tables["bgc_record"]
    edge_params_table = DB.metadata.tables["edge_params"]

    bgc_record_a = alias(bgc_record_table)
    bgc_record_b = alias(bgc_record_table)
    gbk_a = alias(gbk_table)
    gbk_b = alias(gbk_table)

    select_statement = (
        select(
            gbk_a.c.path,
            bgc_record_a.c.record_type,
            bgc_record_a.c.record_number,
            gbk_b.c.path,
            bgc_record_b.c.record_type,
            bgc_record_b.c.record_number,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.ext_a_start,
            distance_table.c.ext_a_stop,
            distance_table.c.ext_b_start,
            distance_table.c.ext_b_stop,
            edge_params_table.c.weights,
            edge_params_table.c.alignment_mode,
            edge_params_table.c.extend_strategy,
        )
        .join(bgc_record_a, distance_table.c.record_a_id == bgc_record_a.c.id)
        .join(bgc_record_b, distance_table.c.record_b_id == bgc_record_b.c.id)
        .join(gbk_a, bgc_record_a.c.gbk_id == gbk_a.c.id)
        .join(gbk_b, bgc_record_b.c.gbk_id == gbk_b.c.id)
        .join(
            edge_params_table,
            distance_table.c.edge_param_id == edge_params_table.c.id,
        )
        .where(distance_table.c.record_a_id.in_(record_ids))
        .where(distance_table.c.record_b_id.in_(record_ids))
        .where(edge_params_table.c.weights.in_(incl_weights))
        .where(edge_params_table.c.alignment_mode == aln_mode)
        .where(edge_params_table.c.extend_strategy == ext_strat)
    )

    if cutoff is not None:
        select_statement = select_statement.where(distance_table.c.distance < cutoff)
    else:
        # still do not include edges with a distance of 1
        select_statement = select_statement.where(distance_table.c.distance < 1)

    edgelist = set(DB.execute(select_statement).fetchall())

    with open(output_path, "w") as network_file:
        header = (
            "GBK_a\tRecord_Type_a\tRecord_Number_a\tORF_coords_a\tGBK_b\t"
            "Record_Type_b\tRecord_Number_b\tORF_coords_b\tdistance\tjaccard\tadjacency\t"
            "dss\tweights\taligmnent_mode\textend_strategy\n"
        )

        network_file.write(header)

        for (
            gbk_path_a,
            record_type_a,
            record_number_a,
            gbk_path_b,
            record_type_b,
            record_number_b,
            distance,
            jaccard,
            adjacency,
            dss,
            ext_a_start,
            ext_a_stop,
            ext_b_start,
            ext_b_stop,
            weights,
            alignment_mode,
            extend_strategy,
        ) in edgelist:
            row = "\t".join(
                [
                    Path(gbk_path_a).stem,
                    record_type_a,
                    str(record_number_a),
                    f"{ext_a_start}:{ext_a_stop}",
                    Path(gbk_path_b).stem,
                    record_type_b,
                    str(record_number_b),
                    f"{ext_b_start}:{ext_b_stop}",
                    f"{distance:.2f}",
                    f"{jaccard:.2f}",
                    f"{adjacency:.2f}",
                    f"{dss:.2f}",
                    weights,
                    alignment_mode,
                    extend_strategy,
                ]
            )

            network_file.write(row + "\n")
