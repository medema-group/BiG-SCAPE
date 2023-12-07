"""Contains config class and method to parse config file """

import configparser
from pathlib import Path


# config class
class BigscapeConfig:
    # static properties

    # PROFILER
    PROFILER_UPDATE_INTERVAL = 0.5

    # LCS
    REGION_MIN_LCS_LEN = 3
    PROTO_MIN_LCS_LEN = 0

    # EXPAND
    REGION_MIN_EXPAND_LEN = 5
    PROTO_MIN_EXPAND_LEN = 3
    NO_MIN_CLASSES = ["Terpene", "Halogenated"]
    EXPAND_MATCH_SCORE = 5
    EXPAND_MISMATCH_SCORE = -3
    EXPAND_GAP_SCORE = -2
    EXPAND_MAX_MATCH_PERC = 0.1

    # TREE
    TOP_FREQS = 3

    @staticmethod
    def parse_config(run) -> None:
        """parses config file

        Args:
            config_file_path (Path): path to config file
        """

        config_file_path = run["config_file_path"]

        config = configparser.ConfigParser()
        config.read(config_file_path)

        # PROFILER
        BigscapeConfig.PROFILER_UPDATE_INTERVAL = float(
            config["PROFILER"]["PROFILER_UPDATE_INTERVAL"]
        )

        # LCS
        BigscapeConfig.REGION_MIN_LCS_LEN = int(config["LCS"]["REGION_MIN_LCS_LEN"])
        BigscapeConfig.PROTO_MIN_LCS_LEN = int(config["LCS"]["PROTO_MIN_LCS_LEN"])

        # EXPAND
        BigscapeConfig.REGION_MIN_EXPAND_LEN = int(
            config["EXPAND"]["REGION_MIN_EXPAND_LEN"]
        )
        BigscapeConfig.PROTO_MIN_EXPAND_LEN = int(
            config["EXPAND"]["PROTO_MIN_EXPAND_LEN"]
        )
        BigscapeConfig.NO_MIN_CLASSES = config["LCS"]["NO_MIN_CLASSES"].split(",")
        BigscapeConfig.EXPAND_MATCH_SCORE = int(config["EXPAND"]["EXPAND_MATCH_SCORE"])
        BigscapeConfig.EXPAND_MISMATCH_SCORE = int(
            config["EXPAND"]["EXPAND_MISMATCH_SCORE"]
        )
        BigscapeConfig.EXPAND_GAP_SCORE = int(config["EXPAND"]["EXPAND_GAP_SCORE"])
        BigscapeConfig.EXPAND_MAX_MATCH_PERC = float(
            config["EXPAND"]["EXPAND_MAX_MATCH_PERC"]
        )

        # TREE
        BigscapeConfig.TOP_FREQS = int(config["TREE"]["TOP_FREQS"])

        # write config log
        BigscapeConfig.write_config_log(run, config)

    @staticmethod
    def write_config_log(run, config) -> None:
        """_summary_

        Args:
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
