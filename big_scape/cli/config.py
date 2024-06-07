"""Contains config class and method to parse config file """

# from python
import yaml
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

        with open(config_file_path) as f:
            config = yaml.load(f, Loader=yaml.FullLoader)

        # PROFILER
        BigscapeConfig.PROFILER_UPDATE_INTERVAL = config["PROFILER_UPDATE_INTERVAL"]

        # INPUT
        BigscapeConfig.MERGED_CAND_CLUSTER_TYPE = config["MERGED_CAND_CLUSTER_TYPE"]
        BigscapeConfig.MIN_BGC_LENGTH = config["MIN_BGC_LENGTH"]
        BigscapeConfig.MAX_BGC_LENGTH = config["MAX_BGC_LENGTH"]

        # LCS
        BigscapeConfig.REGION_MIN_LCS_LEN = config["REGION_MIN_LCS_LEN"]
        BigscapeConfig.PROTO_MIN_LCS_LEN = config["PROTO_MIN_LCS_LEN"]

        # EXPAND
        BigscapeConfig.REGION_MIN_EXPAND_LEN = config["REGION_MIN_EXPAND_LEN"]
        BigscapeConfig.REGION_MIN_EXPAND_LEN_BIO = config["REGION_MIN_EXPAND_LEN_BIO"]
        BigscapeConfig.PROTO_MIN_EXPAND_LEN = config["PROTO_MIN_EXPAND_LEN"]
        BigscapeConfig.NO_MIN_CLASSES = config["NO_MIN_CLASSES"]
        BigscapeConfig.EXPAND_MATCH_SCORE = config["EXPAND_MATCH_SCORE"]
        BigscapeConfig.EXPAND_MISMATCH_SCORE = config["EXPAND_MISMATCH_SCORE"]
        BigscapeConfig.EXPAND_GAP_SCORE = config["EXPAND_GAP_SCORE"]
        BigscapeConfig.EXPAND_MAX_MATCH_PERC = config["EXPAND_MAX_MATCH_PERC"]

        # CLUSTER
        BigscapeConfig.PREFERENCE = config["PREFERENCE"]

        # TREE
        BigscapeConfig.TOP_FREQS = config["TOP_FREQS"]

        # write config log
        BigscapeConfig.write_config_log(run, config)

    @staticmethod
    def write_config_log(run: dict, config: dict) -> None:
        """writes config log file

        Args:
            run (dict): run parameters
            config (configparser.ConfigParser): config settings
        """
        log_path = run["log_path"]
        config_log_path = Path(str(log_path).replace(".log", ".config.log"))

        with open(config_log_path, "w") as config_log:
            for key, value in config.items():
                config_log.write(f"{key}: {value}\n")
