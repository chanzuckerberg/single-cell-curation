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


def retry(func):
    def wrapper(*arg, **kw):
        for i in range(10):
            try:
                res = func(*arg, **kw)
                break
            except Exception:
                print(f"Failed {func.__name__} on try {i+1}")
                time.sleep(2)
        return res

    return wrapper


def log_tracking():
    logging.info("Function: # calls, avg duration, min duration, max duration, total duration")
    for k, v in tracker.items():
        logging.info(f"{k}: {len(v)}, {sum(v)/len(v)}, {min(v)}, {max(v)}, {sum(v)}")
