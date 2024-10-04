"""Contains config class and method to parse config file """

# from python
import yaml
import hashlib
from pathlib import Path
from typing import Optional


# config class
class BigscapeConfig:
    # static default properties
    HASH: str = ""

    # PROFILER
    PROFILER_UPDATE_INTERVAL: float = 0.5

    # INPUT
    MERGED_CAND_CLUSTER_TYPE: list[str] = ["chemical_hybrid", "interleaved"]
    MIN_BGC_LENGTH: int = 0
    MAX_BGC_LENGTH: int = 500000

    # CDS and DOMAIN
    CDS_OVERLAP_CUTOFF: float = 0.1
    DOMAIN_OVERLAP_CUTOFF: float = 0.1

    # LCS
    REGION_MIN_LCS_LEN: float = 0.1
    PROTO_MIN_LCS_LEN: float = 0.0

    # EXTEND
    REGION_MIN_EXTEND_LEN: float = 0.3
    REGION_MIN_EXTEND_LEN_BIO: float = 0.2
    PROTO_MIN_EXTEND_LEN: float = 0.2
    NO_MIN_CLASSES: list[str] = ["terpene"]
    EXTEND_MATCH_SCORE: int = 5
    EXTEND_MISMATCH_SCORE: int = -3
    EXTEND_GAP_SCORE: int = -2
    EXTEND_MAX_MATCH_PERC: float = 0.1

    # CLUSTER
    PREFERENCE: float = 0.0

    # TREE
    TOP_FREQS: int = 3

    # ANCHOR DOMAINS
    ANCHOR_DOMAINS = [
        "PF02801",
        "PF02624",
        "PF00109",
        "PF00501",
        "PF02797",
        "PF01397",
        "PF03936",
        "PF00432",
        "PF00195",
        "PF00494",
        "PF00668",
        "PF05147",
    ]

    # LEGACY ANTISMASH CLASSES
    LEGACY_ANTISMASH_CLASSES = {
        "pks1_products": {"t1pks", "T1PKS"},
        "pksother_products": {
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
        },
        "nrps_products": {"nrps", "NRPS", "NRPS-like", "thioamide-NRP", "NAPAA"},
        "ripps_products": {
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
            "redox-cofactor",
            "thioamitides",
            "epipeptide",
            "cyclic-lactone-autoinducer",
            "spliceotide",
            "RRE-containing",
        },
        "saccharide_products": {
            "amglyccycl",
            "oligosaccharide",
            "cf_saccharide",
            "saccharide",
        },
        "others_products": {
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
        },
    }

    @staticmethod
    def parse_config(config_file_path: Path, log_path: Optional[Path] = None) -> None:
        """parses config file and writes a config.log if log_path is given

        Args:
            config_file_path (Path): path to passed config file
            log_file_path (Optional[Path]): path to log file. Defaults to None.
        """
        print("PARSING CONFIG")
        with open(config_file_path, "rb") as f:
            content = f.read()
            BigscapeConfig.HASH = hashlib.sha256(content).hexdigest()
            config = yaml.load(content, Loader=yaml.FullLoader)

        # PROFILER
        BigscapeConfig.PROFILER_UPDATE_INTERVAL = config["PROFILER_UPDATE_INTERVAL"]

        # INPUT
        BigscapeConfig.MERGED_CAND_CLUSTER_TYPE = config["MERGED_CAND_CLUSTER_TYPE"]
        BigscapeConfig.MIN_BGC_LENGTH = config["MIN_BGC_LENGTH"]
        BigscapeConfig.MAX_BGC_LENGTH = config["MAX_BGC_LENGTH"]

        # CDS and DOMAIN
        BigscapeConfig.CDS_OVERLAP_CUTOFF = config["CDS_OVERLAP_CUTOFF"]
        BigscapeConfig.DOMAIN_OVERLAP_CUTOFF = config["DOMAIN_OVERLAP_CUTOFF"]

        # LCS
        BigscapeConfig.REGION_MIN_LCS_LEN = config["REGION_MIN_LCS_LEN"]
        BigscapeConfig.PROTO_MIN_LCS_LEN = config["PROTO_MIN_LCS_LEN"]

        # EXTEND
        BigscapeConfig.REGION_MIN_EXTEND_LEN = config["REGION_MIN_EXTEND_LEN"]
        BigscapeConfig.REGION_MIN_EXTEND_LEN_BIO = config["REGION_MIN_EXTEND_LEN_BIO"]
        BigscapeConfig.PROTO_MIN_EXTEND_LEN = config["PROTO_MIN_EXTEND_LEN"]
        BigscapeConfig.NO_MIN_CLASSES = config["NO_MIN_CLASSES"]
        BigscapeConfig.EXTEND_MATCH_SCORE = config["EXTEND_MATCH_SCORE"]
        BigscapeConfig.EXTEND_MISMATCH_SCORE = config["EXTEND_MISMATCH_SCORE"]
        BigscapeConfig.EXTEND_GAP_SCORE = config["EXTEND_GAP_SCORE"]
        BigscapeConfig.EXTEND_MAX_MATCH_PERC = config["EXTEND_MAX_MATCH_PERC"]

        # CLUSTER
        BigscapeConfig.PREFERENCE = config["PREFERENCE"]

        # TREE
        BigscapeConfig.TOP_FREQS = config["TOP_FREQS"]

        # ANCHOR DOMAINS
        BigscapeConfig.ANCHOR_DOMAINS = config["ANCHOR_DOMAINS"]

        # LEGACY ANTISMASH CLASSES
        legacy_classes = config["LEGACY_ANTISMASH_CLASSES"]
        for group, classes in legacy_classes.items():
            if isinstance(classes, list):
                legacy_classes[group] = set(classes)
        BigscapeConfig.LEGACY_ANTISMASH_CLASSES = legacy_classes

        # write config log
        if log_path is not None:
            BigscapeConfig.write_config_log(log_path, config)

    @staticmethod
    def write_config_log(log_path: Path, config: dict) -> None:
        """writes config log file

        Args:
            log_path (Path): path to log file
            config (configparser.ConfigParser): config settings
        """
        config_log_path = Path(str(log_path).replace(".log", ".config.log"))

        with open(config_log_path, "w") as config_log:
            for key, value in config.items():
                config_log.write(f"{key}: {value}\n")
