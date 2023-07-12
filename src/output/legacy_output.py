"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
import json
import shutil
from distutils import dir_util
from pathlib import Path
from typing import Any

# from other modules
from src.genbank import GBK, CDS, SOURCE_TYPE
from src.network import BSNetwork


def copy_output_templates(
    output_dir: Path,
    label: str,
    cutoffs: list[float],
    bins: list[str],
):
    """Copy the necessary output html/javascript templates to a relevant output
    directory

    Args:
        output_dir (Path): main output directory
        label (str): label for the run
        cutoffs (list[float]): list of cutoffs used in the analysis
        bins (list[str]): list of bins used in the analyis
    """

    template_root = Path("src/output/html_template")
    template_dir = template_root / "output"
    overview_template = template_root / "overview_html"
    bin_template = template_root / "index_html"

    # copy html content
    dir_util.copy_tree(str(template_dir), str(output_dir))

    # networks subfolders
    output_network_root = output_dir / "html_content/networks"

    if not output_network_root.exists():
        output_network_root.mkdir(exist_ok=True)

    for cutoff in cutoffs:
        cutoff_path = output_network_root / f"{label}_c{cutoff}"

        cutoff_path.mkdir(exist_ok=True)

        # copy overview html
        shutil.copy(str(overview_template), str(cutoff_path / "overview.html"))

        # copy bin html
        for bin in bins:
            bin_path = cutoff_path / bin

            bin_path.mkdir(exist_ok=True)

            shutil.copy(str(bin_template), str(bin_path / "index.html"))


# one of the few direct copy-and-pastes!
def get_class(product):
    """Sort BGC by its type. Uses AntiSMASH annotations
    (see https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)"""

    # TODO: according with current (2021-05) antiSMASH rules:
    # prodigiosin and PpyS-KS -> PKS
    # CDPS -> NRPS
    pks1_products = {"t1pks", "T1PKS"}
    pksother_products = {
        "transatpks",
        "t2pks",
        "t3pks",
        "otherks",
        "hglks",
        "transAT-PKS",
        "transAT-PKS-like",
        "T2PKS",
        "T3PKS",
        "PKS-like",
        "hglE-KS",
    }
    nrps_products = {"nrps", "NRPS", "NRPS-like", "thioamide-NRP", "NAPAA"}
    ripps_products = {
        "lantipeptide",
        "thiopeptide",
        "bacteriocin",
        "linaridin",
        "cyanobactin",
        "glycocin",
        "LAP",
        "lassopeptide",
        "sactipeptide",
        "bottromycin",
        "head_to_tail",
        "microcin",
        "microviridin",
        "proteusin",
        "lanthipeptide",
        "lipolanthine",
        "RaS-RiPP",
        "fungal-RiPP",
        "TfuA-related",
        "guanidinotides",
        "RiPP-like",
        "lanthipeptide-class-i",
        "lanthipeptide-class-ii",
        "lanthipeptide-class-iii",
        "lanthipeptide-class-iv",
        "lanthipeptide-class-v",
        "ranthipeptide",
        "redox-cofactor" "thioamitides" "epipeptide",
        "cyclic-lactone-autoinducer",
        "spliceotide",
        "RRE-containing",
    }
    saccharide_products = {
        "amglyccycl",
        "oligosaccharide",
        "cf_saccharide",
        "saccharide",
    }
    others_products = {
        "acyl_amino_acids",
        "arylpolyene",
        "aminocoumarin",
        "ectoine",
        "butyrolactone",
        "nucleoside",
        "melanin",
        "phosphoglycolipid",
        "phenazine",
        "phosphonate",
        "other",
        "cf_putative",
        "resorcinol",
        "indole",
        "ladderane",
        "PUFA",
        "furan",
        "hserlactone",
        "fused",
        "cf_fatty_acid",
        "siderophore",
        "blactam",
        "fatty_acid",
        "PpyS-KS",
        "CDPS",
        "betalactone",
        "PBDE",
        "tropodithietic-acid",
        "NAGGN",
        "halogenated",
        "pyrrolidine",
    }

    # PKS_Type I
    if product in pks1_products:
        return "PKSI"
    # PKS Other Types
    elif product in pksother_products:
        return "PKSother"
    # NRPs
    elif product in nrps_products:
        return "NRPS"
    # RiPPs
    elif product in ripps_products:
        return "RiPPs"
    # Saccharides
    elif product in saccharide_products:
        return "Saccharides"
    # Terpenes
    elif product == "terpene":
        return "Terpene"
    # PKS/NRP hybrids
    elif len(product.split(".")) > 1:
        # print("  Possible hybrid: (" + cluster + "): " + product)
        # cf_fatty_acid category contains a trailing empty space
        subtypes = set(s.strip() for s in product.split("."))
        if len(subtypes - (pks1_products | pksother_products | nrps_products)) == 0:
            if len(subtypes - nrps_products) == 0:
                return "NRPS"
            elif len(subtypes - (pks1_products | pksother_products)) == 0:
                return "PKSother"  # pks hybrids
            else:
                return "PKS-NRP_Hybrids"
        elif len(subtypes - ripps_products) == 0:
            return "RiPPs"
        elif len(subtypes - saccharide_products) == 0:
            return "Saccharide"
        else:
            return "Others"  # other hybrid
    # Others
    elif product in others_products:
        return "Others"
    # ??
    elif product == "":
        # No product annotation. Perhaps not analyzed by antiSMASH
        return "Others"
    else:
        print("  Warning: unknown product '{}'".format(product))
        return "Others"


def generate_run_data_js(
    output_dir: Path,
    label: str,
    gbks: list[GBK],
    cutoff: float,
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
        gbks (list[GBK]): full list of GBKs used in the analysis
        cutoff (float): cutoff to generate the run_data.js for
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
    families_members: dict[int, list[int]] = {}

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
        region_class = get_class(gbk.region.product)
        if region_class not in class_idx:
            class_idx[region_class] = len(class_idx)

        region_class_idx = class_idx[region_class]

        region_id = str(gbk.path.name)

        # add to list of bgc
        region_bgc_id = len(run_data["input"]["bgc"])
        run_data["input"]["bgc"].append(
            {"acc": region_acc_idx, "class": region_class_idx, "id": region_id}
        )

        if cutoff not in gbk.region._families:
            continue

        family_id = gbk.region._families[cutoff]

        if family_id not in families_members:
            families_members[family_id] = []

        families_members[family_id].append(region_bgc_id)

    # now we have accession list, bgc and families. we can ignore accession_newick for
    # now. construct the classes as the last run_data.input fields
    # order of dict keys guaranteed since python 3.7
    run_data["input"]["classes"] = [
        {"label": classkey} for classkey in list(class_idx.keys())
    ]

    # finally add families
    # TODO: bins
    run_data["networks"].append(
        {
            "label": "mix",
            "families": [],
            "families_similarity": [],
        }
    )
    for family_idx, family_members in families_members.items():
        run_data["networks"][0]["families"].append(
            {
                "label": f"FAM_{family_idx:0>5}",
                "members": family_members,
                "mibig": [],
            }
        )

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


def generate_bigscape_results_js(
    output_dir: Path, label: str, cutoff: float, bins: list[str]
) -> None:
    """Generates the bigscape_results.js found under output/html_content/js

    This function will open an existing file and append a new result, just like v1

    structure of bigscape_results.js:


    {
        "networks": [
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
        bins (list[str]): bins to generate bigscape_results for
    """

    bigscape_results_js_path = output_dir / "html_content/js/bigscape_results.js"

    if bigscape_results_js_path.exists():
        bigscape_results = read_bigscape_results_js(bigscape_results_js_path)
    else:
        bigscape_results = []

    # check if run is already there. if so we can assume its the same run and we won't
    # do anything. this is really only relevant for testing, as runs are usually
    # generated with a timestamp as a label
    cutoff_label = f"{label}_c{cutoff:.1f}"
    for bigscape_result in bigscape_results:
        if bigscape_result["label"] == cutoff_label:
            return

    result_networks = []
    for network in bins:
        if network == "mix":
            result_networks.append(
                {
                    "css": "Others",
                    "label": "Mixed",
                    "name": "mix",
                }
            )
            continue

        result_networks.append(
            {
                "css": network,
                "label": network,
                "name": network,
            }
        )

    bigscape_results.append({"label": cutoff_label, "networks": result_networks})

    with open(bigscape_results_js_path, "w") as run_data_js:
        run_data_js.write(
            "var bigscape_results={};\n".format(
                json.dumps(
                    bigscape_results, indent=4, separators=(",", ":"), sort_keys=True
                )
            )
        )


def generate_bs_data_js_orfs_domains(cds: CDS) -> Any:
    domains = []
    for hsp in cds.hsps:
        domains.append(
            {
                "bitscore": hsp.score,
                "code": hsp.domain,
                "start": hsp.env_start,
                "end": hsp.env_stop,
            }
        )
    return domains


def generate_bs_data_js_orfs(gbk: GBK) -> Any:
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
    gbks: list[GBK],
    cutoff: float,
    bin: str,
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
        gbks (list[GBK]): Full list of GBKs used in the analysis
        cutoff (float): cutoff to generate results for
        bin (str): bin to generate results for
    """
    output_network_root = output_dir / Path("html_content/networks")
    cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")
    bin_path = cutoff_path / bin

    bs_data = []

    # go through all gbks. no need to go through the network
    for gbk in gbks:
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
    gbks: list[GBK], network: BSNetwork
) -> list[list[float]]:
    sim_matrix = []

    for idx, gbk in enumerate(gbks):
        region_sim = []
        region_a = gbk.region
        for other_idx in range(idx):
            region_b = gbks[other_idx].region
            dist = network.graph.adj[region_a][region_b]["dist"]
            sim = round(1 - dist, 4)
            region_sim.append(sim)
        region_sim.append(1.0)  # for self similarity
        sim_matrix.append(region_sim)

    return sim_matrix


def generate_bs_networks_families(
    gbks: list[GBK], cutoff: float
) -> list[dict[str, Any]]:
    # TODO: repeated code can be merged
    networks_families: list[dict[str, Any]] = []
    families_members: dict[int, list[int]] = {}
    for idx, gbk in enumerate(gbks):
        if gbk.region is None:
            raise ValueError("Missing region on GBK!")

        if cutoff not in gbk.region._families:
            continue

        family_id = gbk.region._families[cutoff]

        if family_id not in families_members:
            families_members[family_id] = []

        families_members[family_id].append(idx)

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
    network: BSNetwork,
    gbks: list[GBK],
    cutoff: float,
    bin: str,
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
    bin_path = cutoff_path / bin

    bs_networks_js_path = bin_path / "bs_networks.js"

    # TODO: replace with functions to generate objects
    bs_similarity: list[Any] = generate_bs_networks_js_sim_matrix(gbks, network)
    bs_families: list[Any] = generate_bs_networks_families(gbks, cutoff)
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


def generate_legacy_output(
    output_dir: Path,
    label: str,
    cutoffs: list[float],
    bins: list[str],
    network: BSNetwork,
    gbks: list[GBK],
) -> None:
    copy_output_templates(output_dir, label, cutoffs, bins)

    for cutoff in cutoffs:
        generate_bigscape_results_js(
            output_dir,
            label,
            cutoff,
            bins,
        )
        generate_run_data_js(output_dir, label, gbks, cutoff)

        for bin in bins:
            generate_bs_data_js(output_dir, label, gbks, cutoff, bin)
            generate_bs_networks_js(output_dir, label, network, gbks, cutoff, bin)
