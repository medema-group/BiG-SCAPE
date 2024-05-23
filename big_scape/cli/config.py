"""Contains config class and method to parse config file """

# from python
import configparser
from pathlib import Path


# config class
class BigscapeConfig:
    # static properties

    # PROFILER
    PROFILER_UPDATE_INTERVAL: float = 0.5

    # INPUT
    MERGED_CAND_CLUSTER_TYPE: list[str] = ["chemical_hybrid", "interleaved"]
    MIN_BGC_LENGTH: int = 0
    MAX_BGC_LENGTH: int = 500000

    # LCS
    REGION_MIN_LCS_LEN: int = 3
    PROTO_MIN_LCS_LEN: int = 3

    # EXPAND
    REGION_MIN_EXPAND_LEN: int = 5
    REGION_MIN_EXPAND_LEN_BIO: int = 5
    PROTO_MIN_EXPAND_LEN: int = 3
    NO_MIN_CLASSES: list[str] = ["Terpene"]
    EXPAND_MATCH_SCORE: int = 5
    EXPAND_MISMATCH_SCORE: int = -3
    EXPAND_GAP_SCORE: int = -2
    EXPAND_MAX_MATCH_PERC: float = 0.1

    # CLUSTER
    PREFERENCE: float = 0.0

    # TREE
    TOP_FREQS: int = 3

    @staticmethod
    def parse_config(run: dict) -> None:
        """parses config file

        Args:
            run (dict): run parameters
        """

        config_file_path = run["config_file_path"]

        config = configparser.ConfigParser()
        config.read(config_file_path)

        # PROFILER
        BigscapeConfig.PROFILER_UPDATE_INTERVAL = float(
            config["PROFILER"]["PROFILER_UPDATE_INTERVAL"]
        )

        # INPUT
        BigscapeConfig.MERGED_CAND_CLUSTER_TYPE = config["INPUT"][
            "MERGED_CAND_CLUSTER_TYPE"
        ].split(",")
        BigscapeConfig.MIN_BGC_LENGTH = int(config["INPUT"]["MIN_BGC_LENGTH"])
        BigscapeConfig.MAX_BGC_LENGTH = int(config["INPUT"]["MAX_BGC_LENGTH"])

        # LCS
        BigscapeConfig.REGION_MIN_LCS_LEN = int(config["LCS"]["REGION_MIN_LCS_LEN"])
        BigscapeConfig.PROTO_MIN_LCS_LEN = int(config["LCS"]["PROTO_MIN_LCS_LEN"])

        # EXPAND
        BigscapeConfig.REGION_MIN_EXPAND_LEN = int(
            config["EXPAND"]["REGION_MIN_EXPAND_LEN"]
        )
        BigscapeConfig.REGION_MIN_EXPAND_LEN_BIO = int(
            config["EXPAND"]["REGION_MIN_EXPAND_LEN_BIO"]
        )
        BigscapeConfig.PROTO_MIN_EXPAND_LEN = int(
            config["EXPAND"]["PROTO_MIN_EXPAND_LEN"]
        )
        BigscapeConfig.NO_MIN_CLASSES = config["EXPAND"]["NO_MIN_CLASSES"].split(",")
        BigscapeConfig.EXPAND_MATCH_SCORE = int(config["EXPAND"]["EXPAND_MATCH_SCORE"])
        BigscapeConfig.EXPAND_MISMATCH_SCORE = int(
            config["EXPAND"]["EXPAND_MISMATCH_SCORE"]
        )
        BigscapeConfig.EXPAND_GAP_SCORE = int(config["EXPAND"]["EXPAND_GAP_SCORE"])
        BigscapeConfig.EXPAND_MAX_MATCH_PERC = float(
            config["EXPAND"]["EXPAND_MAX_MATCH_PERC"]
        )

        # CLUSTER
        BigscapeConfig.PREFERENCE = float(config["CLUSTER"]["PREFERENCE"])

        # TREE
        BigscapeConfig.TOP_FREQS = int(config["TREE"]["TOP_FREQS"])

        # write config log
        BigscapeConfig.write_config_log(run, config)

    @staticmethod
    def write_config_log(run: dict, config: configparser.ConfigParser) -> None:
        """writes config log file

        Args:
            run (dict): run parameters
            config (configparser.ConfigParser): config settings
        """
        log_path = run["log_path"]
        config_log_path = Path(str(log_path).replace(".log", ".config.log"))

        with open(config_log_path, "w") as config_log:
            for section in config.sections():
                config_log.write(f"\n[{section}]\n")
                for key in config[section]:
                    value = config[section][key]
                    key = key.upper()
                    config_log.write(f"{key} = {value}\n")
