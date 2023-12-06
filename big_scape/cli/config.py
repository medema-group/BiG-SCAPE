"""Contains config class and method to parse config file """

import configparser


# config class
class BigscapeConfig:
    # static properties

    # PROFILER
    PROFILER_UPDATE_INTERVAL = 0.5

    # LCS
    MIN_LCS_LEN = 3
    NO_MIN_CLASSES = ["Terpene", "Halogenated"]

    # EXPAND
    MIN_EXPAND_LEN = 5
    EXPAND_MATCH_SCORE = 5
    EXPAND_MISMATCH_SCORE = -3
    EXPAND_GAP_SCORE = -2
    EXPAND_MAX_MATCH_PERC = 0.1

    @staticmethod
    def parse_config(config_file_path) -> None:
        config = configparser.ConfigParser()
        config.read(config_file_path)

        # PROFILER
        BigscapeConfig.PROFILER_UPDATE_INTERVAL = float(
            config["PROFILER"]["PROFILER_UPDATE_INTERVAL"]
        )

        # LCS
        BigscapeConfig.MIN_LCS_LEN = int(config["LCS"]["MIN_LCS_LEN"])
        BigscapeConfig.NO_MIN_CLASSES = config["LCS"]["NO_MIN_CLASSES"].split(",")

        # EXPAND
        BigscapeConfig.MIN_EXPAND_LEN = int(config["EXPAND"]["MIN_EXPAND_LEN"])
        BigscapeConfig.EXPAND_MATCH_SCORE = int(config["EXPAND"]["EXPAND_MATCH_SCORE"])
        BigscapeConfig.EXPAND_MISMATCH_SCORE = int(
            config["EXPAND"]["EXPAND_MISMATCH_SCORE"]
        )
        BigscapeConfig.EXPAND_GAP_SCORE = int(config["EXPAND"]["EXPAND_GAP_SCORE"])
        BigscapeConfig.EXPAND_MAX_MATCH_PERC = float(
            config["EXPAND"]["EXPAND_MAX_MATCH_PERC"]
        )
