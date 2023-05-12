import logging
import time
from collections import defaultdict


def configure_logging():
    logging.basicConfig(level=logging.INFO)


tracker = defaultdict(list)


def logit(func):
    def wrapper(*arg, **kw):
        """Logging the start and finish of a function"""
        func_name = func.__name__
        start = time.time()
        # logging.info(f"Start {func_name}", extra={"type": "METRIC"})
        res = func(*arg, **kw)
        duration = time.time() - start
        # logging.info(f"Complete {func_name}", extra={"type": "METRIC"})
        tracker[func_name].append(duration)
        return res

    return wrapper


def print_tracking():
    print("Function: # calls, avg duration, min duration, max duration")
    for k, v in tracker.items():
        print(f"{k}: {len(v)}, {sum(v)/len(v)}, {min(v)}, {max(v)}")
