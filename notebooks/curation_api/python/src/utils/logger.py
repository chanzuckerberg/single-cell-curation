import logging
import os
import sys


def set_log_level(log_level: str):
    """
    :param log_level: the logging level ("NOTSET", "DEBUG", "INFO", "WARN", "ERROR", "FATAL")
    """
    all_levels = logging._nameToLevel.keys()
    if log_level not in all_levels:
        raise Exception(f"The log_level arg must be one of {[x for x in all_levels]}")
    os.environ["log_level"] = log_level
    logger = logging.getLogger()
    logger.setLevel(log_level)
    for h in logger.handlers:
        h.setLevel(log_level)
    print(f"Set logging level to {log_level}")


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
    level_printout = f"{'%(levelname)s:' if logger.level in ('WARN', 'ERROR') else ''}"
    formatter = logging.Formatter(f"{level_printout}%(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    return logger


def failure(logger: logging.Logger, e: Exception, message=None):
    logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m")  # 'FAILED' in bold red
    if message:
        logger.error(message)
    raise e


def success(logger: logging.Logger, message=None):
    logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
    if message:
        logger.info(message)
