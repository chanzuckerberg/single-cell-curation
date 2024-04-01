import json
import logging
import multiprocessing
import random
import tempfile
import threading
import time
from logging.handlers import QueueHandler, QueueListener
from typing import Any, Dict, List
from urllib.request import urlopen

import anndata as ad
import boto3
import numpy as np
import scipy
from anndata.compat import OverloadedDict
from cellxgene_schema.validate import validate
from pandas import DataFrame
from tqdm.contrib.logging import logging_redirect_tqdm

random.seed(42)

MAX_SIZE = 1e9  # Only test datasets that are less than 1GB
SAMPLES = 0.2  # Test 20% of the datasets


def dataset_ids_to_test() -> list[str]:
    logging.info("Fetching dataset ids to test")
    # list all of the datasets
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/curation/v1/datasets").read().decode("utf-8"))
    # # Using a local fiel because the above endpoint is slow
    # with open("./dataset_index.json") as fp:
    #     datasets = json.load(fp)
    logging.info(f"Found {len(datasets)} datasets")
    sample_size = int(len(datasets) * SAMPLES)
    logging.info(f"Sampling {sample_size} datasets")
    filtered_datasets = [dataset for dataset in datasets if dataset["assets"][0]["filesize"] < MAX_SIZE]
    filtered_datasets = sorted(filtered_datasets, key=lambda x: x["assets"][0]["filesize"])
    samples = random.sample(filtered_datasets, sample_size)
    sample_ids = [dataset["dataset_version_id"] for dataset in samples]
    return sample_ids


def iterate_dataset_metadata(dataset_ids) -> List[Dict[str, Any]]:
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/dp/v1/datasets/index").read().decode("utf-8"))
    filtered_dataset_ids = []
    for dataset in datasets:
        if dataset["id"] in dataset_ids:
            filtered_dataset_ids.append(dataset)
    return filtered_dataset_ids


def task(args):
    dataset, log_queue = args
    queue_handler = QueueHandler(log_queue)
    logger = logging.getLogger(f"dataset_{dataset['id']}")
    logger.addHandler(queue_handler)
    logger.setLevel(logging.INFO)

    def download(uri, local_path):
        # Download raw to a temp file
        s3 = boto3.client("s3")
        bucket, key = uri.split("//", 1)[1].split("/", 1)
        s3.download_file(bucket, key, local_path)

    def download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path):
        logger.debug("Downloading")
        raw_h5ad = [asset for asset in dataset["dataset_assets"] if asset["filetype"] == "RAW_H5AD"][0]
        labeled_h5ad = [asset for asset in dataset["dataset_assets"] if asset["filetype"] == "H5AD"][0]

        # Download raw
        t_raw = threading.Thread(target=download, args=(raw_h5ad["s3_uri"], raw_h5ad_path))
        # Download original
        t_labeled = threading.Thread(target=download, args=(labeled_h5ad["s3_uri"], current_labeled_h5ad_path))

        t_raw.start()
        t_labeled.start()
        t_raw.join()
        t_labeled.join()

    errors = []

    def compare_dicts(dict1, dict2, path=None):
        """
        This function recursive compares two dictionaries handling cases where the values
        are unordered arrays with elements that could be dictionaries.
        """
        _errors = []
        path = path if path else []
        dict_types = (dict, OverloadedDict)
        array_types = (np.ndarray, scipy.sparse._csr.csr_matrix)
        if len(dict1) != len(dict2):
            _errors.append(f"{'.'.join(path)} Length of dictionaries not equal")
        for key in dict1:
            path.append(key)
            if key not in dict2:
                _errors.append(f"{'.'.join(path)} Key {key} not in dict2")
                continue
            _match = True
            value1 = dict1[key]
            value2 = dict2[key]

            if isinstance(value1, dict_types) and isinstance(value2, dict_types):
                if not compare_dicts(value1, value2, path):
                    _errors.append(f"{'.'.join(path)} {key} does not match")
            elif isinstance(value1, list) and isinstance(value2, list):
                if len(value1) != len(value2):
                    _errors.append(f"{'.'.join(path)} List lengths do not match")
                # check if the lists contain dictionaries as elements
                if len(value1) > 0 and isinstance(value1[0], dict_types) and isinstance(value2[0], dict_types):
                    for i in range(len(value1)):
                        if not compare_dicts(value1[i], value2[i], path):
                            _errors.append(f"{'.'.join(path)} Lists do not match")
                            break
                elif sorted(value1) != sorted(value2):
                    _errors.append(f"{'.'.join(path)} Lists do not match")
            elif isinstance(value1, array_types) and type(value1) == type(value2):
                if not np.array_equal(value1, value2):
                    _errors.append(f"{'.'.join(path)} arrays do not match ")
            elif isinstance(value1, DataFrame) and type(value1) == type(value2):
                if not value1.equals(value2):
                    _errors.append(f"{'.'.join(path)} DataFrames do not match")
            else:
                try:
                    if value1 != value2:
                        _errors.append(f"{'.'.join(path)} values do not match")
                except ValueError:
                    logger.exception(type(value1))
            path.pop()
        if _errors:
            errors.extend(_errors)
            return False
        return True

    try:
        with tempfile.TemporaryDirectory() as path:
            current_labeled_h5ad_path = f"{path}/{dataset['id']}_current_labeled.h5ad"
            raw_h5ad_path = f"{path}/{dataset['id']}_raw.h5ad"
            new_labeled_h5ad_path = f"{path}/{dataset['id']}_new_labeled.h5ad"

            download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path)
            logger.debug("Validating")
            validate(raw_h5ad_path, new_labeled_h5ad_path)

            # Compare h5ads
            logger.debug("Comparing")
            new_ad = ad.read_h5ad(new_labeled_h5ad_path, backed="r")
            current_ad = ad.read_h5ad(current_labeled_h5ad_path, backed="r")
            if not current_ad.obs.equals(new_ad.obs):
                errors.append("obs not equal")
            if not current_ad.var.equals(new_ad.var):
                errors.append("var not equal")
            new_ad.uns["citation"] = current_ad.uns["citation"]
            compare_dicts(current_ad.uns, new_ad.uns)
            if errors:
                logger.error("FAILED\n\t" + "\n\t".join(errors))
            else:
                logger.info("PASS")
    except Exception as e:
        logger.exception(f"Error testing dataset: {e}")


def main():
    # configure logging
    mp_manager = multiprocessing.Manager()
    log_queue = mp_manager.Queue()
    queue_handler = QueueHandler(log_queue)

    logger = logging.getLogger()
    logger.addHandler(queue_handler)
    logger.setLevel(logging.INFO)

    formatter = logging.Formatter("%(asctime)s:%(process)d:%(name)s:%(levelname)s:%(message)s")

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(f"test_cog_with_dataset_index_{int(time.time())}.log")
    file_handler.setFormatter(formatter)

    listener = QueueListener(log_queue, console_handler, file_handler)
    listener.start()

    # get dataset ids to test
    test_datasets_ids = dataset_ids_to_test()
    dataset_ids = "\n\t".join(test_datasets_ids)
    logger.info(f"Datasets Versions:\n\t{dataset_ids}")

    with logging_redirect_tqdm():
        # from tqdm import tqdm
        # for dataset in tqdm(
        #         iterate_dataset_metadata(test_dataset_ids)):
        #     task([dataset, log_queue])
        from tqdm.contrib import concurrent

        concurrent.process_map(
            task, [(d, log_queue) for d in iterate_dataset_metadata(test_datasets_ids)], max_workers=4
        )

    listener.stop()


if __name__ == "__main__":
    main()
