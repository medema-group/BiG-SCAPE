"""Contains functions to mimic legacy output as seen in BiG-SCAPE 1.0"""

# from python
import json
import shutil
from distutils import dir_util
from pathlib import Path
from typing import Any, Optional

# from other modules
from src.genbank import GBK
from src.network import BSNetwork


def copy_output_templates(
    output_dir: Path, label: str, cutoffs: Optional[list[float]] = None
):
    """Copy the necessary output html/javascript templates to a relevant output
    directory
    """

    if cutoffs is None:
        cutoffs = []

    template_root = Path("src/output/html_template")
    template_dir = template_root / Path("output")
    overview_template = template_root / Path("overview_html")

    # copy html content
    dir_util.copy_tree(str(template_dir), str(output_dir))

    # networks subfolders
    output_network_root = output_dir / Path("html_content/networks")

    if not output_network_root.exists():
        output_network_root.mkdir(exist_ok=True)

    for cutoff in cutoffs:
        cutoff_path = output_network_root / Path(f"{label}_c{cutoff}")

        cutoff_path.mkdir(exist_ok=True)

        # copy overview html
        shutil.copy(str(overview_template), str(cutoff_path))


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


def generate_run_data(
    output_dir: Path, network: BSNetwork, gbks: list[GBK], cutoff: float, label: str
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
    """Reads an existing bigscape_results_js into a dictionary"""

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
    output_dir: Path, label: str, cutoff: float, networks: list[str]
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
    for network in networks:
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


def generate_legacy_output(
    output_dir: Path,
    label: str,
    cutoffs: list[float],
    bins: list[str],
    network: BSNetwork,
    gbks: list[GBK],
) -> None:
    copy_output_templates(output_dir, label, cutoffs)

    for cutoff in cutoffs:
        generate_bigscape_results_js(output_dir, label, cutoff, bins)
        generate_run_data(output_dir, network, gbks, cutoff, label)
