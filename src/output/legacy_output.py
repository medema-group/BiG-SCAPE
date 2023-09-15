"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
import json
import shutil
from distutils import dir_util
from pathlib import Path
from typing import Any

# from other modules
from src.comparison import RecordPairGenerator, legacy_get_class
from src.genbank import GBK, CDS
from src.network import BSNetwork
from src.enums import SOURCE_TYPE


def copy_base_output_templates(output_dir: Path):
    """Copy the base output html/javascript templates to an output
    directory

    Args:
        output_dir (Path): main output directory
    """

    template_root = Path("src/output/html_template")
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

    template_root = Path("src/output/html_template")
    overview_template = template_root / "overview_html"

    # copy overview html
    shutil.copy(str(overview_template), str(cutoff_path / "overview.html"))


def prepare_bin_folder(
    output_dir: Path, label: str, cutoff: float, bin: RecordPairGenerator
) -> None:
    """Prepare the output folder for a bin under a cutoff

    Args:
        output_dir (Path): output folder
        label (str): run label
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin
    """
    # networks subfolders
    output_network_root = output_dir / "html_content/networks"

    cutoff_path = output_network_root / f"{label}_c{cutoff}"

    bin_path = cutoff_path / bin.label

    bin_path.mkdir(exist_ok=True)

    template_root = Path("src/output/html_template")
    bin_template = template_root / "index_html"

    shutil.copy(str(bin_template), str(bin_path / "index.html"))


def generate_pfams_js(output_dir: Path, pfam_info: list[tuple[str, str, str]]) -> None:
    """Generate the pfam.js file needed to show the correct PFAM domain colors

    Args:
        output_dir (Path): main output directory
        pfam_info (list[tuple[str]]): A list of tuples containing pfam information. Each tuple
        contains an accession, a name and a description
    """

    # gather color information
    pfam_colors_file_path = Path("src/output/domain_colors.tsv")

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
        color = "255,255,255"
        if accession in pfam_colors_dict:
            color = pfam_colors_dict[accession]

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
    output_dir: Path,
    label: str,
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
            "label": string, (e.g. PKSother. these are the bins. can also be "mix")
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

    run_data: dict[str, Any] = {
        "duration": "Test",
        "start_time": "Test",
        "end_time": "Test",
        "parameters": "Test",
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
        region_class = legacy_get_class(gbk.region.product)
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
    bin: RecordPairGenerator,
    families_members: dict[int, list[int]],
) -> None:
    """Add run data to the run_data.js file for the given bin with a given cutoff

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff
        bin (BGCBin): BGC bin
        families_members (dict[int, list[int]]): a dictionary with family ids as keys
        and a list of region indexes as values
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    run_data_js_path = cutoff_path / Path("run_data.js")

    run_data = read_run_data_js(run_data_js_path)

    run_data["networks"].append(
        {
            "label": bin.label,
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
    Networks are appended to later when each bin is generated

    structure of bigscape_results.js:


    {
        "networks": [ // per bin
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
    output_dir: Path, label: str, cutoff: float, bin: RecordPairGenerator
) -> None:
    """Add cutoff and network data to an existing bigscape_results.js file

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin

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
    css = bin.label
    label = bin.label

    # one exception to the above values is the mix bin, which uses others as a CSS class
    if bin.label == "mix":
        css = "Others"
        label = "Mixed"

    bigscape_result["networks"].append(
        {
            "css": css,
            "label": label,
            "name": bin.label,
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
    bin: RecordPairGenerator,
):
    """Generates the bs_data.js file located at
    output/html_content/networks/[label]/[bin]

    structure:
    [
        {
            "desc": str, (e.g. Streptomyces coelicolor A3(2) complete genome)
            "start: int,
            "end": int,
            "id": str, (e.g. AL645882.2.cluster010),
            "mibig": bool,
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
        bin (str): bin to generate results for
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    bin_path = cutoff_path / bin.label

    bs_data = []

    # go through all gbks. no need to go through the network
    for record in bin.source_records:
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
                "orfs": generate_bs_data_js_orfs(gbk),
            }
        )

    bs_data_path = bin_path / "bs_data.js"

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
    bin: RecordPairGenerator,
    network: BSNetwork,
) -> list[list[float]]:
    """Generate a similarity matrix for the bs_networks.js file

    Args:
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin
        network (BSNetwork): network object

    Returns:
        list[list[float]]: a similarity matrix
    """
    sim_matrix = []

    for idx, record_a in enumerate(bin.source_records):
        region_sim = []
        for record_b in bin.source_records[:idx]:
            dist = network.graph.adj[record_a][record_b]["dist"]
            if dist > cutoff:
                dist = 1.0
            sim = round(1 - dist, 4)
            region_sim.append(sim)
        region_sim.append(1.0)  # for self similarity
        sim_matrix.append(region_sim)

    return sim_matrix


def generate_bs_families_members(
    cutoff: float,
    bin: RecordPairGenerator,
) -> dict[int, list[int]]:
    """Generate a dictionary where keys are family indexes and values are list of region
    indexes that belong to that family

    Args:
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin

    Returns:
        dict[int, list[int]]: family to member regions index
    """
    families_members: dict[int, list[int]] = {}
    for idx, record in enumerate(bin.source_records):
        if cutoff not in record._families:
            continue

        family_id = record._families[cutoff]

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

    Returns:
        list[dict[str, Any]]: networks families object
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
    bin: RecordPairGenerator,
    network: BSNetwork,
    bs_families: list[dict[str, Any]],
):
    """Generates the bs_networks.js file located at
    output/html_content/networks/[label]/[bin]

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
        bin (str): bin to generate results for
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    bin_path = cutoff_path / bin.label

    bs_networks_js_path = bin_path / "bs_networks.js"

    # TODO: replace with functions to generate objects
    bs_similarity: list[Any] = generate_bs_networks_js_sim_matrix(cutoff, bin, network)
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
    to bins

    Args:
        output_dir (Path): output directory
        pfam_info (list[tuple(str, str, str)]): A list of tuples containing information
        about the PFAM entries used in this analysis. tuples correspond to the accession
        name and description of a pfam entry
    """
    copy_base_output_templates(output_dir)

    generate_pfams_js(output_dir, pfam_info)


def legacy_prepare_cutoff_output(
    output_dir: Path, label: str, cutoff: float, gbks: list[GBK]
) -> None:
    """Prepare output data for a given cutoff value

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        gbks (list[GBK]): list of gbks used in the analysis
    """
    prepare_cutoff_folder(output_dir, label, cutoff)

    generate_bigscape_results_js(output_dir, label, cutoff)

    generate_run_data_js(output_dir, label, cutoff, gbks)


def legacy_prepare_bin_output(
    output_dir: Path, label: str, cutoff: float, bin: RecordPairGenerator
) -> None:
    """Prepare output data for a given bin at a given cutoff value

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin
    """
    prepare_bin_folder(output_dir, label, cutoff, bin)
    generate_bs_data_js(output_dir, label, cutoff, bin)
    add_bigscape_results_js_network(output_dir, label, cutoff, bin)


def legacy_generate_bin_output(
    output_dir: Path,
    label: str,
    cutoff: float,
    bin: RecordPairGenerator,
    network: BSNetwork,
) -> None:
    """Generate the network data from a bin from cutoff filtering and affinity
    propagation

    Args:
        output_dir (Path): output directory
        label (str): run label
        cutoff (float): cutoff value
        bin (BGCBin): BGC bin
        network (BSNetwork): the network object for the bin
    """
    families_members = generate_bs_families_members(cutoff, bin)
    networks_families = generate_bs_networks_families(families_members)

    add_run_data_network(output_dir, label, cutoff, bin, families_members)
    generate_bs_networks_js(output_dir, label, cutoff, bin, network, networks_families)