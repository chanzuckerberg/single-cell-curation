import logging
import sys
import typing


def get_custom_logger(log_level: typing.Union[str, int]) -> logging.Logger:
    """
    Get a custom logger that will still print to stdout in notebooks.
    :param log_level: the logging level ("INFO" or 20, "WARN" or 30, etc.)
    :return: the logger object
    """
    logging.basicConfig(level=log_level)
    logger = logging.getLogger()
    logger.removeHandler(logger.handlers[0])
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(level=log_level)
    formatter = logging.Formatter("%(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger
