import logging
import os
import sys


def set_log_level(log_level: str):
    """
    :param log_level: the logging level ("NOTSET", "DEBUG", "INFO", "WARN", "ERROR", "FATAL")
    """
    all_levels = logging._nameToLevel.keys()
    if log_level not in all_levels:
        raise Exception(f"The log_level arg must be one of {list(all_levels)}")
    os.environ["LOG_LEVEL"] = log_level
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
    log_level = os.getenv("LOG_LEVEL", "INFO")
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


def failure(logger: logging.Logger, e, *messages):
    logger.error("\n\033[1m\033[38;5;9mFAILED\033[0m\n")  # 'FAILED' in bold red
    logger.error(e.response.reason + "\n")
    logger.error(e.response.text + "\n")
    logger.error("x-request-id: {}".format(e.response.headers.get("x-request-id")))

    if messages:
        [logger.error(m) for m in messages]


def success(logger: logging.Logger, *messages):
    logger.info("\n\033[1m\033[38;5;10mSUCCESS\033[0m\n")  # 'SUCCESS' in bold green
    if messages:
        [logger.info(m) for m in messages]
