import logging
import os
import sys


def set_log_level(log_level: str = "INFO"):
    """
    :param log_level: the logging level ("INFO" or 20, "WARN" or 30, etc...)
    """
    os.environ["log_level"] = log_level


def get_custom_logger() -> logging.Logger:
    """
    Get a custom logger that will still print to stdout in notebooks.
    :return: the logger object
    """
    log_level = os.getenv("log_level", "INFO")
    logging.basicConfig(level=log_level)
    logger = logging.getLogger()
    logger.removeHandler(logger.handlers[0])
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level=log_level)
    level_printout = f"{'%(levelname)s:' if log_level in ('WARN', 'ERROR') else ''}"
    formatter = logging.Formatter(f"{level_printout}%(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger
