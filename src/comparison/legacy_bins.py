"""Contains code to generate legacy bins based on product."""


# one of the few direct copy-and-pastes!
def legacy_get_class(product):
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
