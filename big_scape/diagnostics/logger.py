"""Module containing code to handle logging"""

# from python
import logging
import sys


def init_logger(run) -> None:  # pragma: no cover
    """Initialize the logger"""

    # logger
    # this tells the logger what the messages should look like
    # asctime = YYYY-MM-DD HH:MM:SS,fff
    # levelname = DEBUG/INFO/WARN/ERROR
    # message = whatever we pass, eg logging.debug("message")
    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")

    # get the built in logger
    root_logger = logging.getLogger()

    if run["verbose"]:
        root_logger.level = logging.DEBUG
    else:
        root_logger.level = logging.INFO

    if not run["quiet"]:
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(log_formatter)
        root_logger.addHandler(console_handler)


def init_logger_file(run) -> None:  # pragma: no cover
    """Initialize the logger file"""

    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")
    root_logger = logging.getLogger()
    file_handler = logging.FileHandler(run["log_path"])
    file_handler.setFormatter(log_formatter)
    root_logger.addHandler(file_handler)
    logging.info(" ".join(sys.argv))
