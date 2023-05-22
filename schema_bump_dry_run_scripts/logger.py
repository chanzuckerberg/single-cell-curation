import logging
from collections import defaultdict


def configure_logging():
    logging.basicConfig(level=logging.INFO)


tracker = defaultdict(list)
