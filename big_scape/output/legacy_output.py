"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
import json
import shutil
from distutils import dir_util
from pathlib import Path
from typing import Any
from sqlalchemy import select

# from other modules
from big_scape.data import DB
from big_scape.comparison import RecordPairGenerator, legacy_get_class
from big_scape.genbank import GBK, CDS, BGCRecord
from big_scape.enums import SOURCE_TYPE
from big_scape.comparison import get_record_category

import big_scape.network.network as bs_network
import big_scape.network.utility as bs_network_utility


def copy_base_output_templates(output_dir: Path):
    """Copy the base output html/javascript templates to an output
    directory

    Args:
        output_dir (Path): main output directory
    """

    template_root = Path("big_scape/output/html_template")
    template_dir = template_root / "output"

    # copy html content
    dir_util.copy_tree(str(template_dir), str(output_dir), verbose=False)

    # networks subfolders
    output_network_root = output_dir / "html_content/networks"

    if not output_network_root.exists():
        output_network_root.mkdir(exist_ok=True)


def prepare_cutoff_folder(output_dir: Path, label: str, cutoff: float) -> None:
    """Prepare a folder for a given cutoff output

    Args:
        output_dir (Path): base output directory
        label (str): run label
        cutoff (float): cutoff value
    """
    # networks subfolders
    output_network_root = output_dir / "html_content/networks"

    cutoff_path = output_network_root / f"{label}_c{cutoff}"

    cutoff_path.mkdir(exist_ok=True)

    template_root = Path("big_scape/output/html_template")
    overview_template = template_root / "overview_html"

    # copy overview html
    shutil.copy(str(overview_template), str(cutoff_path / "overview.html"))


def prepare_pair_generator_folder(
    output_dir: Path, label: str, cutoff: float, pair_generator: RecordPairGenerator
) -> None:
    """Prepare the output folder for a pair_generator under a cutoff

    Args:
        output_dir (Path): output folder
        label (str): run label
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator
    """
    # networks subfolders
    output_network_root = output_dir / "html_content/networks"

    cutoff_path = output_network_root / f"{label}_c{cutoff}"

    pair_generator_path = cutoff_path / pair_generator.label

    pair_generator_path.mkdir(exist_ok=True)

    template_root = Path("big_scape/output/html_template")
    pair_generator_template = template_root / "index_html"

    shutil.copy(str(pair_generator_template), str(pair_generator_path / "index.html"))


def generate_pfams_js(output_dir: Path, pfam_info: list[tuple[str, str, str]]) -> None:
    """Generate the pfam.js file needed to show the correct PFAM domain colors

    Args:
        output_dir (Path): main output directory
        pfam_info (list[tuple[str]]): A list of tuples containing pfam information. Each tuple
        contains an accession, a name and a description
    """

    # gather color information
    pfam_colors_file_path = Path("big_scape/output/domain_colors.tsv")

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


def generate_run_data_js(
    run: dict,
    cutoff: float,
    gbks: list[GBK],
):
    """Generate the run_data.js file needed for the html output

    structure:

    "duration": string,
    "start_time": string.
    "end_time": string,
    "parameters": string,
    "input": {
        "accession": [
            {
                "id": string,
                "label": string,
            }
        ],
        "accession_newick": ?[],
        "bgc": [
            {
                "acc": int,
                "class": int, (refers to one of the class_idx below)
                "id": string,
            }
        ]
    },
    "networks": [
        {
            "label": string, (e.g. PKSother. these are the pair_generators. can also be "mix")
            "families": [
                "label": "FAM_xxxxx",
                "members": int[], (refers to bgc index above),
                "mibig": ?[]
            ],
            "families_similarity": [[[]]], (some sort of similarity matrix),
        }
    ]

    Args:
        output_dir (Path): main output path
        label (str): label of the run
        cutoff (float): cutoff to generate the run_data.js for
        gbks (list[GBK]): full list of GBKs used in the analysis
    """

    output_dir = run["output_dir"]
    label = run["label"]

    run_data: dict[str, Any] = {
        "duration": str(run["duration"])[:-7],
        "start_time": str(run["start_time"])[:-7],
        "end_time": str(run["end_time"])[:-7],
        "mode": run["mode"],
        "input_dir": str(run["input_dir"].name) if run["input_dir"] else "None",
        "output_dir": str(run["output_dir"].name) if run["output_dir"] else "None",
        "reference_dir": str(run["reference_dir"].name)
        if run["reference_dir"]
        else "None",
        "query_path": str(run["query_bgc_path"].name)
        if run["query_bgc_path"]
        else "NA",
        "mibig": run["mibig_version"] if run["mibig_version"] else "None",
        "record_type": run["record_type"].name.title(),
        "min_bgc_length": run["min_bgc_length"],
        "classify": "Legacy Groups"
        if run["legacy_classify"]
        else run["classify"].name.title()
        if run["classify"]
        else "Not Classify",
        "weights": "Legacy Weights" if run["legacy_weights"] else "Mix",
        "alignment_mode": run["alignment_mode"].name.title(),
        "include_singletons": "NA"
        if "include_singletons" not in run.keys()
        else ("Yes" if run["include_singletons"] else "No"),
        "input": {
            "accession": [],
            "accession_newick": [],
            "bgc": [],
        },
        "networks": [],
    }

    # these are mostly index dictionaries needed for certain fields
    members: dict[GBK, int] = {}
    genomes: dict[str, int] = {}
    class_idx: dict[str, int] = {}

    for idx, gbk in enumerate(gbks):
        if gbk.region is None:
            continue

        # add to idx dict
        members[gbk] = len(members)

        # build accs
        if "organism" in gbk.metadata:
            organism = gbk.metadata["organism"]
        else:
            organism = "Unknown"

        if organism not in genomes:
            # set to new id of run data accessions
            genomes[organism] = len(run_data["input"]["accession"])
            # add to run dat accessions
            run_data["input"]["accession"].append(
                {
                    "id": f"genome_{genomes[organism]}",
                    "label": organism,
                }
            )

        # acc id of gbk
        region_acc_idx = genomes[organism]

        # class id
        # product hybrids of AS4 and under dealt with here and in legacy_bin_generator
        product = ".".join(gbk.region.product.split("-"))
        region_class = legacy_get_class(product)

        if region_class not in class_idx:
            class_idx[region_class] = len(class_idx)

        region_class_idx = class_idx[region_class]

        region_id = str(gbk.path.name)

        # add to list of bgc
        run_data["input"]["bgc"].append(
            {"acc": region_acc_idx, "class": region_class_idx, "id": region_id}
        )

        if cutoff not in gbk.region._families:
            continue

    # now we have accession list, bgc and families. we can ignore accession_newick for
    # now. construct the classes as the last run_data.input fields
    # order of dict keys guaranteed since python 3.7
    run_data["input"]["classes"] = [
        {"label": classkey} for classkey in list(class_idx.keys())
    ]

    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    run_data_js_path = cutoff_path / Path("run_data.js")

    with open(run_data_js_path, "w") as run_data_js:
        run_data_js.write(
            "var run_data={};\n".format(
                json.dumps(run_data, indent=4, separators=(",", ":"), sort_keys=True)
            )
        )
        run_data_js.write("dataLoaded();\n")


def add_run_data_network(
    output_dir: Path,
    label: str,
    cutoff: float,
    pair_generator: RecordPairGenerator,
    families_members: dict[int, list[int]],
) -> None:
    """Add run data to the run_data.js file for the given pair_generator with a given cutoff

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff
        pair_generator (BGCBin): BGC pair_generator
        families_members (dict[int, list[int]]): a dictionary with family ids as keys
        and a list of region indexes as values
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    run_data_js_path = cutoff_path / Path("run_data.js")

    run_data = read_run_data_js(run_data_js_path)

    run_data["networks"].append(
        {
            "label": pair_generator.label,
            "families": [],
            "families_similarity": [],
        }
    )

    for family_idx, family_members in families_members.items():
        run_data["networks"][-1]["families"].append(
            {
                "label": f"FAM_{family_idx:0>5}",
                "members": family_members,
                "mibig": [],
            }
        )

    with open(run_data_js_path, "w") as run_data_js:
        run_data_js.write(
            "var run_data={};\n".format(
                json.dumps(run_data, indent=4, separators=(",", ":"), sort_keys=True)
            )
        )
        run_data_js.write("dataLoaded();\n")


def read_bigscape_results_js(bigscape_results_js_path: Path) -> list[Any]:
    """Reads an existing bigscape_results_js into a dictionary by stripping the
    JavaScript parts and decoding the JSON content

    Args:
        bigscape_results_js_path (Path): path to the existing bigscape_results.js

    Returns:
        list[Any]: the bigscape_results.js content as an object
    """

    lines = []
    with open(bigscape_results_js_path, mode="r", encoding="utf-8") as bigscape_results:
        # first line has the "var bigscape_Results = " bit we want to get rid of
        first_line = bigscape_results.readline()
        bracket_idx = first_line.index("[")
        remove_left = bracket_idx
        lines.append(first_line[remove_left:])

        lines.extend(bigscape_results.readlines())

    # last one has a semicolon at the end. remove it
    lines[-1] = lines[-1][:1]

    return json.loads("".join(lines))


def read_run_data_js(run_data_js_path: Path) -> dict[str, Any]:
    """Reads an existing run_data.js into a dictionary by stripping the
    JavaScript parts and decoding the JSON content

    Args:
        run_data_js_path (Path): path to the existing run_data.js

    Returns:
        dict[Any]: the run_data.js content as an object
    """

    lines = []
    with open(run_data_js_path, mode="r", encoding="utf-8") as bigscape_results:
        # first line has the "var bigscape_Results = " bit we want to get rid of
        first_line = bigscape_results.readline()
        bracket_idx = first_line.index("{")
        remove_left = bracket_idx
        lines.append(first_line[remove_left:])

        lines.extend(bigscape_results.readlines())

    # last line is the dataLoaded(); statement. remove it
    lines.pop()

    # last line after that has a semicolon at the end. remove the semicolon
    lines[-1] = lines[-1][:1]

    return json.loads("".join(lines))


def generate_bigscape_results_js(output_dir: Path, label: str, cutoff: float) -> None:
    """Generates the bigscape_results.js found under output/html_content/js

    This function will open an existing file and append a new result, just like v1
    If a result with a label already exists, the networks in that run are cleared.
    Networks are appended to later when each pair_generator is generated

    structure of bigscape_results.js:


    {
        "networks": [ // per pair_generator
            {
                "css": string, (e.g. PKSother. mix uses 'Others'),
                "label": string, (appears at the top of the overview as a tab),
                "name": string, (appears at the network overview as a tab)
            }
        ]
    }

    Args:
        output_dir (Path): main output directory
        label (str): label of the run
        cutoff (float): cutoff to generate bigscape_results for
    """

    bigscape_results_js_path = output_dir / "html_content/js/bigscape_results.js"

    if bigscape_results_js_path.exists():
        bigscape_results = read_bigscape_results_js(bigscape_results_js_path)
    else:
        bigscape_results = []

    # check if run is already there. if so we can re use it
    cutoff_label = f"{label}_c{cutoff:.1f}"
    bigscape_result = None
    for existing_result in bigscape_results:
        if existing_result["label"] == cutoff_label:
            bigscape_result = existing_result
            # we will reset the networks list in case it has changed
            bigscape_result["networks"] = []
            break

    if bigscape_result is None:
        bigscape_results.append({"label": cutoff_label, "networks": []})

        bigscape_result = bigscape_results[-1]

    with open(bigscape_results_js_path, "w") as run_data_js:
        run_data_js.write(
            "var bigscape_results={};\n".format(
                json.dumps(
                    bigscape_results, indent=4, separators=(",", ":"), sort_keys=True
                )
            )
        )


def add_bigscape_results_js_network(
    output_dir: Path, label: str, cutoff: float, pair_generator: RecordPairGenerator
) -> None:
    """Add cutoff and network data to an existing bigscape_results.js file

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator

    Raises:
        FileNotFoundError: Raised when the bigscape_results.js file does not exist
    """
    bigscape_results_js_path = output_dir / "html_content/js/bigscape_results.js"

    if not bigscape_results_js_path.exists():
        raise FileNotFoundError("bigscape_results.js not found when it should exist!")

    bigscape_results = read_bigscape_results_js(bigscape_results_js_path)

    # check if run is already there. if so we can re use it
    cutoff_label = f"{label}_c{cutoff:.1f}"
    bigscape_result: dict[str, Any] = {"label": cutoff_label, "networks": []}

    for existing_result in bigscape_results:
        if existing_result["label"] == cutoff_label:
            bigscape_result = existing_result
            break

    # the following is a bit confusing because the mix class uses some different values
    css = pair_generator.label
    label = pair_generator.label

    # one exception to the above values is the mix pair_generator, which uses others as a CSS class
    if pair_generator.label == "mix":
        css = "Others"
        label = "mix"

    bigscape_result["networks"].append(
        {
            "css": css,
            "label": label,
            "name": pair_generator.label,
        }
    )

    with open(bigscape_results_js_path, "w") as run_data_js:
        run_data_js.write(
            "var bigscape_results={};\n".format(
                json.dumps(
                    bigscape_results, indent=4, separators=(",", ":"), sort_keys=True
                )
            )
        )


def generate_bs_data_js_orfs_domains(cds: CDS) -> list[dict[str, Any]]:
    """Generate the domains for the orfs field in a bs_data.js file

    Args:
        cds (CDS): Coding Domain Sequence object

    Returns:
        list[dict[str: Any]]: A list of domains in the form of dictionaries with
        bitscores, accessions, start coordinates and stop coordinates
    """
    domains = []
    for hsp in cds.hsps:
        domains.append(
            {
                "bitscore": hsp.score,
                "code": hsp.domain[:7],  # trim version number
                "start": hsp.env_start,
                "end": hsp.env_stop,
            }
        )
    return domains


def generate_bs_data_js_orfs(gbk: GBK) -> Any:
    """Generate the orf fields for an existing bs_data.js file

    Args:
        gbk (GBK): GBK file

    Returns:
        Any: A dictionary object representing an ORF as understood by the big-scape
        output page
    """
    orfs = []
    for idx, cds in enumerate(gbk.genes):
        orfs.append(
            {
                "id": f"{gbk.path.name[:-4]}_ORF{idx+1}",
                "start": cds.nt_start,
                "end": cds.nt_stop,
                "strand": 1 if cds.strand else -1,
                "domains": generate_bs_data_js_orfs_domains(cds),
            }
        )
    return orfs


def generate_bs_data_js(
    output_dir: Path,
    label: str,
    cutoff: float,
    pair_generator: RecordPairGenerator,
):
    """Generates the bs_data.js file located at
    output/html_content/networks/[label]/[pair_generator]

    structure:
    [
        {
            "desc": str, (e.g. Streptomyces coelicolor A3(2) complete genome)
            "start: int,
            "end": int,
            "id": str, (e.g. AL645882.2.cluster010),
            "mibig": bool,
            "source": str, (e.g. mibig, reference, or query),
            "orfs": [
                {
                    "domains": [
                        {
                            "bitscore": float,
                            "code": str, (e.g. PF08241.14),
                            "start": int, (env start I am guessing)
                            "end": int
                        }
                    ],
                    "id": str, (e.g. AL645882.2.cluster010_ORF1),
                    "start": int,
                    "end": int,
                    "strand": int (1 or -1)
                }
            ]
        }
    ]

    Args:
        output_dir (Path): main output directory
        label (str): _description_
        cutoff (float): cutoff to generate results for
        pair_generator (str): pair_generator to generate results for
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    pair_generator_path = cutoff_path / pair_generator.label

    bs_data = []

    # go through all gbks. no need to go through the network
    for record in pair_generator.source_records:
        if record.parent_gbk is None:
            raise AttributeError("Record parent GBK is not set!")

        gbk = record.parent_gbk
        organism = "Unknown"
        if "organism" in gbk.metadata:
            organism = gbk.metadata["organism"]

        bs_data.append(
            {
                "desc": organism,
                "start": 1,
                "end": len(gbk.nt_seq),
                "id": gbk.path.name,
                "mibig": gbk.source_type == SOURCE_TYPE.MIBIG,
                "source": gbk.source_type.name.lower(),
                "orfs": generate_bs_data_js_orfs(gbk),
            }
        )

    bs_data_path = pair_generator_path / "bs_data.js"

    with open(bs_data_path, "w", encoding="utf-8") as bs_data_js:
        bs_data_js.write(
            "var bs_data={};\ndataLoaded('bs_data');\n".format(
                json.dumps(
                    bs_data,
                    indent=4,
                    separators=(",", ":"),
                    sort_keys=True,
                )
            )
        )


def generate_bs_networks_js_sim_matrix(
    cutoff: float,
    pair_generator: RecordPairGenerator,
) -> list[list[float]]:
    """Generate a similarity matrix for the bs_networks.js file

    Args:
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator
        network (BSNetwork): network object

    Returns:
        list[list[float]]: a similarity matrix
    """
    edges = bs_network.get_edges(pair_generator.record_ids, cutoff)

    adj_list = bs_network_utility.edge_list_to_adj_list(edges)

    # the above adjacency list is likely not ordered in the same way as the entries in
    # the pair_generator, which are eventually used to determine the order of records in the output
    # files. Go through the records in the same order as the pair_generator and create a sparse
    # matrix from that
    # the adjacency list contains edges, also. these contain distances. convert to sim
    # by subracting distance from 1

    sparse_matrix = []

    for record_idx, record_a in enumerate(pair_generator.source_records):
        a_record_id = record_a._db_id

        sparse_row = []
        # the sparse matrix is the 'lower half' of a full similarity matrix.
        # this means that if we are on the third element, matrix so far should look like
        # 1
        # 1 2
        # 1 2 3
        # obviously this means that the last element in each row should always be 1.0
        # and we are iterating "behind" the current element
        for record_b in pair_generator.source_records[:record_idx]:
            # assume similarity is 0.0
            similarity = 0.0

            # but if we find an entry in the adjacency list this can change
            if a_record_id in adj_list:
                if record_b._db_id in adj_list[a_record_id]:
                    similarity = adj_list[a_record_id][record_b._db_id]

            sparse_row.append(similarity)

        # add the 1.0 to the end
        sparse_row.append(1.0)

        sparse_matrix.append(sparse_row)

    return sparse_matrix


def generate_bs_families_members(
    cutoff: float,
    pair_generator: RecordPairGenerator,
) -> dict[int, list[int]]:
    """Generate a dictionary where keys are family indexes and values are list of region
    indexes that belong to that family

    Args:
        cutoff (float): cutoff value
        pair_generator (BGCBin): BGC pair_generator

    Returns:
        dict[int, list[int]]: family to member regions index
    """
    # get a dictionary of node id to family id
    node_family = {}

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    # get all families from the database
    bgc_families_table = DB.metadata.tables["bgc_record_family"]

    select_statement = (
        bgc_families_table.select()
        .where(bgc_families_table.c.record_id.in_(pair_generator.record_ids))
        .where(bgc_families_table.c.cutoff == cutoff)
    )

    result = DB.execute(select_statement).fetchall()

    for row in result:
        node_family[row.record_id] = row.family

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


def generate_bs_networks_families(
    families_members: dict[int, list[int]]
) -> list[dict[str, Any]]:
    """Generate a list object that represents the families that are present in a network
    and the regions that are members of that family for the bs_networks.js file

    Args:
        families_members (dict[int, list[int]]): family to member regions index
    """
    networks_families: list[dict[str, Any]] = []
    for family_idx, family_members in families_members.items():
        networks_families.append(
            {
                "id": f"FAM_{family_idx:0>5}",
                "members": family_members,
            }
        )
    return networks_families


def generate_bs_networks_js(
    output_dir: Path,
    label: str,
    cutoff: float,
    pair_generator: RecordPairGenerator,
    bs_families: list[dict[str, Any]],
):
    """Generates the bs_networks.js file located at
    output/html_content/networks/[label]/[pair_generator]

    There are four components to this file:
    bs_similarity
    bs_families
    bs_families_alignment
    bs_similarity_families

    each have their own function to generate the objects, and are written to a file here

    Since this is annoying to do, we're not going to attempt re-reading and appending
    to an existing file. we just overwrite.

    Args:
        output_dir (Path): main output directory
        label (str): label of the run
        network (BSNetwork): BiG-SCAPE network object of nodes (bgcs/regions) and edges
        gbks (list[GBK]): Full list of GBK files TODO: may not be needed. remove?
        cutoff (float): cutoff to generate results for
        pair_generator (str): pair_generator to generate results for
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    pair_generator_path = cutoff_path / pair_generator.label

    bs_networks_js_path = pair_generator_path / "bs_networks.js"

    # TODO: replace with functions to generate objects
    bs_similarity: list[Any] = generate_bs_networks_js_sim_matrix(
        cutoff, pair_generator
    )
    bs_families_alignment: list[Any] = []
    bs_similarity_families: list[Any] = []

    with open(bs_networks_js_path, "w") as bs_networks_js:
        bs_networks_js.write(
            "var bs_similarity={};\n".format(
                json.dumps(
                    bs_similarity,
                    indent=4,
                    separators=(",", ":"),
                    sort_keys=True,
                )
            )
        )

        bs_networks_js.write(
            "var bs_families={};\n".format(
                json.dumps(
                    bs_families,
                    indent=4,
                    separators=(",", ":"),
                    sort_keys=True,
                )
            )
        )

        bs_networks_js.write(
            "var bs_families_alignment={};\n".format(
                json.dumps(
                    bs_families_alignment,
                    indent=4,
                    separators=(",", ":"),
                    sort_keys=True,
                )
            )
        )

        bs_networks_js.write(
            "var bs_similarity_families={};\n".format(
                json.dumps(
                    bs_similarity_families,
                    indent=4,
                    separators=(",", ":"),
                    sort_keys=True,
                )
            )
        )

        bs_networks_js.write("dataLoaded('bs_networks');\n")


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
    copy_base_output_templates(output_dir)

    generate_pfams_js(output_dir, pfam_info)


def legacy_prepare_cutoff_output(run: dict, cutoff: float, gbks: list[GBK]) -> None:
    """Prepare output data for a given cutoff value

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        gbks (list[GBK]): list of gbks used in the analysis
    """
    prepare_cutoff_folder(run["output_dir"], run["label"], cutoff)

    generate_bigscape_results_js(run["output_dir"], run["label"], cutoff)

    generate_run_data_js(run, cutoff, gbks)


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

    output_dir = run["output_dir"]
    label = run["label"]

    prepare_pair_generator_folder(output_dir, label, cutoff, pair_generator)
    generate_bs_data_js(output_dir, label, cutoff, pair_generator)
    add_bigscape_results_js_network(output_dir, label, cutoff, pair_generator)


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
    output_dir = run["output_dir"]
    label = run["label"]

    families_members = generate_bs_families_members(cutoff, pair_generator)
    networks_families = generate_bs_networks_families(families_members)

    add_run_data_network(output_dir, label, cutoff, pair_generator, families_members)
    generate_bs_networks_js(
        output_dir, label, cutoff, pair_generator, networks_families
    )
    write_clustering_file(run, cutoff, pair_generator)
    write_cutoff_network_file(run, cutoff, pair_generator)


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

    output_dir = run["output_dir"]
    label = run["label"]
    output_network_root = output_dir / "html_content/networks"
    cutoff_path = output_network_root / f"{label}_c{cutoff}"
    record_annotations_path = cutoff_path / "record_annotations.tsv"

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")
    bgc_record_table = DB.metadata.tables["bgc_record"]
    gbk_table = DB.metadata.tables["gbk"]

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

        for record in all_bgc_records:
            record_id = record._db_id

            select_statment = select(
                bgc_record_table.c.gbk_id,
                bgc_record_table.c.record_number,
                bgc_record_table.c.record_type,
                bgc_record_table.c.product,
            ).where(bgc_record_table.c.id == record_id)
            gbk_id, record_number, record_type, product = DB.execute(
                select_statment
            ).fetchone()

            category = get_record_category(record)

            select_statment = select(
                gbk_table.c.path,
                gbk_table.c.organism,
                gbk_table.c.taxonomy,
                gbk_table.c.description,
            ).where(gbk_table.c.id == gbk_id)

            gbk_path, organism, taxonomy, description = DB.execute(
                select_statment
            ).fetchone()
            # gbk_path = gbk_path.split("/")[-1]

            row = "\t".join(
                [
                    str(gbk_path),
                    record_type,
                    str(record_number),
                    product,
                    category,
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

    output_network_root = output_dir / "html_content/networks"
    cutoff_path = output_network_root / f"{label}_c{cutoff}"
    pair_generator_path = cutoff_path / pair_generator.label
    clustering_file_path = pair_generator_path / f"{bin_label}_clustering_c{cutoff}.tsv"

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    gbk_table = DB.metadata.tables["gbk"]
    bgc_record_table = DB.metadata.tables["bgc_record"]
    record_famly_table = DB.metadata.tables["bgc_record_family"]

    record_ids = pair_generator.record_ids

    with open(clustering_file_path, "w") as clustering_file:
        header = "\t".join(
            [
                "GBK",
                "Record_Type",
                "Record_Number",
                "GCF_number",
            ]
        )
        clustering_file.write(header + "\n")

        select_statement = (
            select(
                record_famly_table.c.record_id,
                record_famly_table.c.family,
                record_famly_table.c.cutoff,
            )
            .where(record_famly_table.c.record_id.in_(record_ids))
            .where(record_famly_table.c.cutoff == cutoff)
        )

        rows = DB.execute(select_statement).fetchall()

        for row in rows:
            record_id, gcf_number, cutoff = row

            select_statment = select(
                bgc_record_table.c.gbk_id,
                bgc_record_table.c.record_number,
                bgc_record_table.c.record_type,
            ).where(bgc_record_table.c.id == record_id)
            gbk_id, record_number, record_type = DB.execute(select_statment).fetchone()

            select_statment = select(gbk_table.c.path).where(gbk_table.c.id == gbk_id)
            gbk_path = DB.execute(select_statment).fetchone()[0]

            row = "\t".join(
                [
                    str(gbk_path),
                    record_type,
                    str(record_number),
                    str(gcf_number),
                ]
            )
            clustering_file.write(row + "\n")

    return None


def write_network_file(
    network_file_path: Path,
    existing_distances: set[tuple[int, int, float, float, float, float, str]],
) -> None:
    """Writes a network file to the output directory

    Args:
        network_file_path (Path): path to write file to
        existing_distances (set): edgelist to populate file

    Raises:
        RuntimeError: no db present

    """

    if not DB.metadata:
        raise RuntimeError("DB.metadata is None")

    gbk_table = DB.metadata.tables["gbk"]
    bgc_record_table = DB.metadata.tables["bgc_record"]

    with open(network_file_path, "w") as network_file:
        header = "\t".join(
            [
                "GBK_a",
                "Record_Type_a",
                "Record_Number_a",
                "GBK_b",
                "Record_Type_b",
                "Record_Number_b",
                "distance",
                "jaccard",
                "adjacency",
                "dss",
                "weights",
            ]
        )

        network_file.write(header + "\n")

        for dist in existing_distances:
            region_a_id, region_b_id, distance, jaccard, adjacency, dss, weights = dist

            # for record A
            select_statment = select(
                bgc_record_table.c.gbk_id,
                bgc_record_table.c.record_number,
                bgc_record_table.c.record_type,
            ).where(bgc_record_table.c.id == region_a_id)
            gbk_id_a, record_number_a, record_type_a = DB.execute(
                select_statment
            ).fetchone()

            select_statment = select(gbk_table.c.path).where(gbk_table.c.id == gbk_id_a)
            gbk_path_a = DB.execute(select_statment).fetchone()[0]

            # for record B
            select_statment = select(
                bgc_record_table.c.gbk_id,
                bgc_record_table.c.record_number,
                bgc_record_table.c.record_type,
            ).where(bgc_record_table.c.id == region_b_id)
            gbk_id_b, record_number_b, record_type_b = DB.execute(
                select_statment
            ).fetchone()

            select_statment = select(gbk_table.c.path).where(gbk_table.c.id == gbk_id_b)
            gbk_path_b = DB.execute(select_statment).fetchone()[0]

            row = "\t".join(
                [
                    str(gbk_path_a),
                    record_type_a,
                    str(record_number_a),
                    str(gbk_path_b),
                    record_type_b,
                    str(record_number_b),
                    f"{distance:.2f}",
                    f"{jaccard:.2f}",
                    f"{adjacency:.2f}",
                    f"{dss:.2f}",
                    weights,
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

    output_network_root = output_dir / "html_content/networks"
    cutoff_path = output_network_root / f"{label}_c{cutoff}"
    pair_generator_path = cutoff_path / pair_generator.label
    cutoff_network_file_path = pair_generator_path / f"{bin_label}_c{cutoff}.network"

    # get distances for this bin
    existing_distances = get_cutoff_edgelist(run, cutoff, pair_generator)

    # write file
    write_network_file(cutoff_network_file_path, existing_distances)

    return None


def get_cutoff_edgelist(
    run: dict, cutoff: float, pair_generator: RecordPairGenerator
) -> set[tuple[int, int, float, float, float, float, str]]:
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

    select_statement = (
        select(
            distance_table.c.region_a_id,
            distance_table.c.region_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.weights,
        )
        .where(distance_table.c.region_a_id.in_(pair_generator.record_ids))
        .where(distance_table.c.region_b_id.in_(pair_generator.record_ids))
        .where(distance_table.c.weights == pair_generator.weights)
        .where(distance_table.c.distance < cutoff)
    )

    existing_distances = set(DB.execute(select_statement).fetchall())

    return existing_distances


def write_full_network_file(run: dict, all_bgc_records: list[BGCRecord]) -> None:
    """Writes the full network file to the output directory,
    i.e. all edges from db that have both records in the run,
    and for all weights relevant to the run

    Args:
        run (dict): run parameters
        all_bgc_records (list): BGC records in this run
    """

    output_dir = run["output_dir"]
    label = run["label"]

    output_network_root = output_dir / "html_content/networks"
    full_network_file_path = output_network_root / f"{label}_full.network"

    # get distances for this set of records
    existing_distances = get_full_network_edgelist(run, all_bgc_records)

    # write file
    write_network_file(full_network_file_path, existing_distances)


def get_full_network_edgelist(
    run: dict, all_bgc_records: list
) -> set[tuple[int, int, float, float, float, float, str]]:
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

    select_statement = (
        select(
            distance_table.c.region_a_id,
            distance_table.c.region_b_id,
            distance_table.c.distance,
            distance_table.c.jaccard,
            distance_table.c.adjacency,
            distance_table.c.dss,
            distance_table.c.weights,
        )
        .where(distance_table.c.region_a_id.in_(record_ids))
        .where(distance_table.c.region_b_id.in_(record_ids))
        .where(distance_table.c.weights.in_(incl_weights))
    )

    existing_distances = set(DB.execute(select_statement).fetchall())

    return existing_distances
