"""Module containing code to handle logging"""

# from python
import logging

# from other modules
from src.parameters import RunParameters


def init_logger(run: RunParameters) -> None:
    """Initializes the logger"""

    # logger
    # this tells the logger what the messages should look like
    # asctime = YYYY-MM-DD HH:MM:SS,fff
    # levelname = DEBUG/INFO/WARN/ERROR
    # message = whatever we pass, eg logging.debug("message")
    log_formatter = logging.Formatter("%(asctime)s %(levelname)-7.7s %(message)s")

    # get the built in logger
    root_logger = logging.getLogger()

    if run.diagnostics.verbose:
        root_logger.level = logging.DEBUG
    else:
        root_logger.level = logging.INFO

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(log_formatter)
    root_logger.addHandler(console_handler)
