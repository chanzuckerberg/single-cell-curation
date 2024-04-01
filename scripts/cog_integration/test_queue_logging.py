import logging
import multiprocessing
import time
from logging.handlers import QueueHandler, QueueListener
from multiprocessing import get_context

from tqdm.contrib.concurrent import process_map
from tqdm.contrib.logging import logging_redirect_tqdm


def task(log_queue):
    queue_handler = QueueHandler(log_queue)
    logger = logging.getLogger()
    logger.addHandler(queue_handler)
    logger.setLevel(logging.INFO)

    time.sleep(1)
    logger.info("Task")


def main():
    mp_manager = multiprocessing.Manager()
    log_queue = mp_manager.Queue()
    queue_handler = QueueHandler(log_queue)

    logger = logging.getLogger()
    logger.addHandler(queue_handler)
    logger.setLevel(logging.INFO)

    console_handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s:%(process)d:%(name)s:%(levelname)s:%(message)s")
    console_handler.setFormatter(formatter)

    file_handler = logging.FileHandler("test.log")
    file_handler.setFormatter(formatter)

    listener = QueueListener(log_queue, console_handler, file_handler)

    listener.start()
    with logging_redirect_tqdm():
        logger.info("Starting")
        process_map(task, [log_queue for _ in range(10)], max_workers=4, context=get_context("fork"))
        logger.info("Done")

    listener.stop()


if __name__ == "__main__":
    main()
