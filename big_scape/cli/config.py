"""Contains config class and method to parse config file """

# from python
import configparser
from pathlib import Path


# config class
class BigscapeConfig:
    # static properties

    # PROFILER
    PROFILER_UPDATE_INTERVAL = 0.5

    # INPUT
    MERGED_CAND_CLUSTER_TYPE = ["chemical_hybrid", "interleaved"]
    MIN_BGC_LENGTH = 0
    MAX_BGC_LENGTH = 500000

    # LCS
    REGION_MIN_LCS_LEN = 3
    PROTO_MIN_LCS_LEN = 3

    # EXPAND
    REGION_MIN_EXPAND_LEN = 5
    REGION_MIN_EXPAND_LEN_BIO = 5
    PROTO_MIN_EXPAND_LEN = 3
    NO_MIN_CLASSES = ["Terpene"]
    EXPAND_MATCH_SCORE = 5
    EXPAND_MISMATCH_SCORE = -3
    EXPAND_GAP_SCORE = -2
    EXPAND_MAX_MATCH_PERC = 0.1

    # CLUSTER
    # TODO: benchmark these values, or at least think better about them
    EDGE_WEIGHT_STD_THRESHOLD = 0.1
    CC_CONNECTIVITY_THRESHOLD = 0.8
    BETWEENNESS_CENTRALITY_NODES = 0.3

    # TREE
    TOP_FREQS = 3

    @staticmethod
    def parse_config(run) -> None:
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
        BigscapeConfig.EDGE_WEIGHT_STD_THRESHOLD = float(
            config["CLUSTER"]["EDGE_WEIGHT_STD_THRESHOLD"]
        )
        BigscapeConfig.CC_CONNECTIVITY_THRESHOLD = float(
            config["CLUSTER"]["CC_CONNECTIVITY_THRESHOLD"]
        )
        BigscapeConfig.BETWEENNESS_CENTRALITY_NODES = float(
            config["CLUSTER"]["BETWEENNESS_CENTRALITY_NODES"]
        )

        # TREE
        BigscapeConfig.TOP_FREQS = int(config["TREE"]["TOP_FREQS"])

        # write config log
        BigscapeConfig.write_config_log(run, config)

    @staticmethod
    def write_config_log(run, config) -> None:
        """writes config lof file

        Args:
            run (dict): run parameters
            log_path (Path): path to log file
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
