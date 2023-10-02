"""Contains code to generate legacy bins based on product."""

# from python
from typing import Iterator

# from other modules
from big_scape.genbank import GBK, BGCRecord

# from this module
from .binning import RecordPairGenerator


# weights are in the order JC, AI, DSS, Anchor boost
LEGACY_BINS = {
    "PKSI": {"weights": (0.22, 0.02, 0.76, 1.0)},
    "PKSother": {"weights": (0.0, 0.68, 0.32, 4.0)},
    "NRPS": {"weights": (0.0, 0.0, 1.0, 4.0)},
    "RiPPs": {"weights": (0.28, 0.01, 0.71, 1.0)},
    "Saccharides": {"weights": (0.0, 1.0, 0.0, 1.0)},
    "Terpene": {"weights": (0.2, 0.05, 0.75, 2.0)},
    "PKS-NRP_Hybrids": {"weights": (0.0, 0.22, 0.78, 1.0)},
    "Others": {"weights": (0.01, 0.02, 0.97, 4.0)},
    "mix": {"weights": (0.2, 0.05, 0.75, 2.0)},
}


def legacy_bin_generator(
    gbks: list[GBK],
) -> Iterator[RecordPairGenerator]:  # pragma no cover
    """Generate bins for each class as they existed in the BiG-SCAPE 1.0 implementation

    Args:
        gbks (list[GBK]): List of GBKs to generate bins for

    Yields:
        Iterator[BGCBin]: Generator that yields bins. Order is not guarenteed to be
        consistent
        TODO: should it be consistent?
    """
    # generate index
    class_idx: dict[str, list[BGCRecord]] = {
        class_name: [] for class_name in LEGACY_BINS.keys() if class_name != "mix"
    }

    for gbk in gbks:
        if gbk.region is None:
            continue

        if gbk.region.product is None:
            continue

        region_class = legacy_get_class(gbk.region.product)

        class_idx[region_class].append(gbk.region)

    for class_name, regions in class_idx.items():
        bin = RecordPairGenerator(class_name)
        bin.add_records(regions)
        yield bin


# one of the few direct copy-and-pastes!
def legacy_get_class(product):  # pragma no cover
    """Sort BGC by its type. Uses AntiSMASH annotations
    (see https://docs.antismash.secondarymetabolites.org/glossary/#cluster-types)

    Args:
        product (str): product type

    Returns:
        str: product class
    """

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
