"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
from itertools import repeat
import json
from multiprocessing import Pool
import platform
from distutils import dir_util
from pathlib import Path
import click
from sqlalchemy import select, alias

# from other modules
from big_scape.data import DB
from big_scape.comparison import RecordPairGenerator
from big_scape.genbank import GBK, BGCRecord
from big_scape.trees import generate_newick_tree, save_tree, save_trees
from big_scape.comparison import get_record_category

import big_scape.paths as bs_paths


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

    if click_context and click_context.obj["no_interactive"]:
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

    if click_context and click_context.obj["no_interactive"]:
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

    if click_context and click_context.obj["no_interactive"]:
        return

    output_dir: Path = run["output_dir"]
    label = run["label"]

    # networks subfolders
    output_files_root = output_dir / "output_files"
    cutoff_files_path = output_files_root / f"{label}_c{cutoff}"
    pair_generator_files_path = cutoff_files_path / pair_generator.label
    pair_generator_files_path.mkdir(exist_ok=True)

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

    if click_context and click_context.obj["no_interactive"]:
        return

    families_members = generate_bs_families_members(
        cutoff, pair_generator, run["run_id"]
    )
    write_clustering_file(run, cutoff, pair_generator)
    write_cutoff_network_file(run, cutoff, pair_generator)
    generate_newick_trees(run, cutoff, pair_generator, families_members)


def generate_newick_trees(
    run: dict,
    cutoff: float,
    pair_generator: RecordPairGenerator,
    families_members: dict[int, list[int]],
):
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

    system = platform.system()

    exemplars_db = []
    families_records = []
    exemplars_idx = []
    for exemplar, members in families_members.items():
        family_recs = [pair_generator.source_records[member] for member in members]
        for idx, record in enumerate(family_recs):
            if record._db_id == exemplar:
                exemplar_idx = idx
                break
        if system == "Darwin":
            tree = generate_newick_tree(family_recs, exemplar_idx, exemplar, tree_path)
            save_tree(tree, exemplar, run["run_id"])
        else:
            exemplars_db.append(exemplar)
            families_records.append(family_recs)
            exemplars_idx.append(exemplar_idx)

    if system != "Darwin":
        pool = Pool(processes=run["cores"])
        trees = pool.starmap(
            generate_newick_tree,
            zip(families_records, exemplars_idx, exemplars_db, repeat(tree_path)),
        )
        save_trees(trees, exemplars_db, run["run_id"])
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

    if click_context and click_context.obj["no_interactive"]:
        return

    output_dir = run["output_dir"]
    label = run["label"]
    output_files_root = output_dir / "output_files"
    cutoff_path = output_files_root / f"{label}_c{cutoff}"
    record_annotations_path = cutoff_path / "record_annotations.tsv"

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")
    bgc_record_table = DB.metadata.tables["bgc_record"]
    gbk_table = DB.metadata.tables["gbk"]

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
        .where(bgc_record_table.c.id.in_(record_categories.keys()))
        .join(bgc_record_table, bgc_record_table.c.gbk_id == gbk_table.c.id)
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

    record_ids = pair_generator.record_ids
    select_statement = (
        select(
            gbk_table.c.path,
            bgc_record_table.c.record_type,
            bgc_record_table.c.record_number,
            rec_fam_table.c.family_id,
            family_table.c.center_id,
        )
        .join(bgc_record_table, bgc_record_table.c.gbk_id == gbk_table.c.id)
        .join(rec_fam_table, bgc_record_table.c.id == rec_fam_table.c.record_id)
        .join(family_table, rec_fam_table.c.family_id == family_table.c.id)
        .where(rec_fam_table.c.record_id.in_(record_ids))
        .where(family_table.c.cutoff == cutoff)
        .where(family_table.c.bin_label == bin_label)
    )

    record_data = DB.execute(select_statement).fetchall()

    with open(clustering_file_path, "w") as clustering_file:
        header = "\t".join(
            ["GBK", "Record_Type", "Record_Number", "GCF_Number", "Family_Name"]
        )
        clustering_file.write(header + "\n")

        for record in record_data:
            gbk_path, record_type, record_number, gcf_number, center = record

            row = "\t".join(
                [
                    Path(gbk_path).stem,
                    record_type,
                    str(record_number),
                    str(gcf_number),
                    f"FAM_{center:0>5}",
                ]
            )
            clustering_file.write(row + "\n")

    return None


def write_network_file(
    network_file_path: Path,
    existing_distances: set[
        tuple[
            str,
            str,
            int,
            str,
            str,
            int,
            float,
            float,
            float,
            float,
            str,
            int,
            int,
            int,
            int,
            str,
        ]
    ],
) -> None:
    """Writes a network file to the output directory

    Args:
        network_file_path (Path): path to write file to
        existing_distances (set): edgelist to populate file

    Raises:
        RuntimeError: no db present

    """

    with open(network_file_path, "w") as network_file:
        header = "\t".join(
            [
                "GBK_a",
                "Record_Type_a",
                "Record_Number_a",
                "ORF_coords_a",
                "GBK_b",
                "Record_Type_b",
                "Record_Number_b",
                "ORF_coords_b",
                "distance",
                "jaccard",
                "adjacency",
                "dss",
                "weights",
                "aligmnent_mode",
            ]
        )

        network_file.write(header + "\n")

        # TODO: optimize this, get all info from db in one query without passing existing distances as intermediate step
        for dist in existing_distances:
            (
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
                weights,
                ext_a_start,
                ext_a_stop,
                ext_b_start,
                ext_b_stop,
                alignment_mode,
            ) = dist

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
                ]
            )

            network_file.write(row + "\n")

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

    # get distances for this bin
    existing_distances = get_cutoff_edgelist(run, cutoff, pair_generator)

    # write file
    write_network_file(cutoff_network_file_path, existing_distances)

    return None


def get_cutoff_edgelist(
    run: dict, cutoff: float, pair_generator: RecordPairGenerator
) -> set[
    tuple[
        str,
        str,
        int,
        str,
        str,
        int,
        float,
        float,
        float,
        float,
        str,
        int,
        int,
        int,
        int,
        str,
    ]
]:
    """Generate the network egdelist for a given bin with edges above the cutoff

    Args:
        run (dict): run parameters

    Raises:
        RuntimeError: no database present

    Returns:
        set: edgelist
    """

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
            edge_params_table.c.weights,
            distance_table.c.ext_a_start,
            distance_table.c.ext_a_stop,
            distance_table.c.ext_b_start,
            distance_table.c.ext_b_stop,
            edge_params_table.c.alignment_mode,
        )
        .join(bgc_record_a, distance_table.c.record_a_id == bgc_record_a.c.id)
        .join(bgc_record_b, distance_table.c.record_b_id == bgc_record_b.c.id)
        .join(gbk_a, bgc_record_a.c.gbk_id == gbk_a.c.id)
        .join(gbk_b, bgc_record_b.c.gbk_id == gbk_b.c.id)
        .join(
            edge_params_table,
            distance_table.c.edge_param_id == edge_params_table.c.id,
        )
        .where(distance_table.c.record_a_id.in_(pair_generator.record_ids))
        .where(distance_table.c.record_b_id.in_(pair_generator.record_ids))
        .where(edge_params_table.c.weights == pair_generator.weights)
        .where(distance_table.c.distance < cutoff)
    )

    edgelist = set(DB.execute(select_statement).fetchall())

    return edgelist


def write_full_network_file(run: dict, all_bgc_records: list[BGCRecord]) -> None:
    """Writes the full network file to the output directory,
    i.e. all edges from db that have both records in the run,
    and for all weights relevant to the run

    Args:
        run (dict): run parameters
        all_bgc_records (list): BGC records in this run
    """

    click_context = click.get_current_context(silent=True)

    if click_context and click_context.obj["no_interactive"]:
        return

    output_dir = run["output_dir"]
    label = run["label"]

    output_files_root = output_dir / "output_files"
    full_network_file_path = output_files_root / f"{label}_full.network"

    # get distances for this set of records
    edgelist = get_full_network_edgelist(run, all_bgc_records)

    # write file
    write_network_file(full_network_file_path, edgelist)


def get_full_network_edgelist(
    run: dict, all_bgc_records: list
) -> set[
    tuple[
        str,
        str,
        int,
        str,
        str,
        int,
        float,
        float,
        float,
        float,
        str,
        int,
        int,
        int,
        int,
        str,
    ]
]:
    """Get all edges for the pairs of records in this run, for the weights relevant to this run

    Args:
        run (dict): run parameters

    Raises:
        RuntimeError: no database present

    Returns:
        set: edgelist
    """

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

    if run["no_mix"]:
        incl_weights.remove("mix")
    if run["legacy_weights"]:
        incl_weights.extend(legacy_weights)
    # an empty list should never happen since no_mix and no classify are mutually exclusive in the CLI validations

    record_ids = [record._db_id for record in all_bgc_records]

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
            edge_params_table.c.weights,
            distance_table.c.ext_a_start,
            distance_table.c.ext_a_stop,
            distance_table.c.ext_b_start,
            distance_table.c.ext_b_stop,
            edge_params_table.c.alignment_mode,
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
    )

    edgelist = set(DB.execute(select_statement).fetchall())

    return edgelist
