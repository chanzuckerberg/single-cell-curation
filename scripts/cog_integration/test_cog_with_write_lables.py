import json
import logging
import random
import tempfile
import threading
import time
from typing import Any, Dict, List, Optional
from urllib.request import urlopen

import anndata as ad
import boto3
import numpy as np
from anndata.compat import OverloadedDict
from cellxgene_schema.validate import validate
from tqdm.contrib.logging import logging_redirect_tqdm

# Configure logging
console_log = logging.StreamHandler()
console_log.setLevel(logging.ERROR)
file_log = logging.FileHandler(f"test_cog_with_dataset_index_{int(time.time())}.log")
file_log.setLevel(logging.INFO)
logging.basicConfig(level=logging.ERROR, handlers=[console_log, file_log])
main_logger = logging.getLogger(__name__)

random.seed(42)

MAX_SIZE = 1e9  # Only test datasets that are less than 1GB
SAMPLES = 0.2  # Test 20% of the datasets


def dataset_ids_to_test() -> list[str]:
    main_logger.info("Fetching dataset ids to test")
    # list all of the datasets
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/curation/v1/datasets").read().decode("utf-8"))
    main_logger.info(f"Found {len(datasets)} datasets")
    sample_size = int(len(datasets) * SAMPLES)
    main_logger.info(f"Sampling {sample_size} datasets")
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


def download(uri, local_path, logger):
    # Download raw to a temp file
    s3 = boto3.client("s3")
    logger.info(f"Downloading {uri} to {local_path}")
    bucket, key = uri.split("//", 1)[1].split("/", 1)
    s3.download_file(bucket, key, local_path)
    logger.info(f"Downloaded {uri} to {local_path}")


def download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path, logger):
    raw_h5ad = [asset for asset in dataset["dataset_assets"] if asset["filetype"] == "RAW_H5AD"][0]
    labeled_h5ad = [asset for asset in dataset["dataset_assets"] if asset["filetype"] == "H5AD"][0]

    # Download raw
    t_raw = threading.Thread(target=download, args=(raw_h5ad["s3_uri"], raw_h5ad_path, logger))
    # Download original
    t_labeled = threading.Thread(target=download, args=(labeled_h5ad["s3_uri"], current_labeled_h5ad_path, logger))

    t_raw.start()
    t_labeled.start()
    t_raw.join()
    t_labeled.join()


def compare_dicts(dict1, dict2, logger):
    """
    This function recursive compares two dictionaries handling cases where the values
    are unordered arrays with elements that could be dictionaries.
    """
    dict_types = (dict, OverloadedDict)
    if len(dict1) != len(dict2):
        logger.error("Length of dictionaries not equal")
    match = True
    for key in dict1:
        if key not in dict2:
            logger.error(f"Key {key} not in dict2")
            continue
        _match = True
        value1 = dict1[key]
        value2 = dict2[key]

        if isinstance(value1, dict_types) and isinstance(value2, dict_types):
            if not compare_dicts(value1, value2, logger):
                _match = False
        elif isinstance(value1, list) and isinstance(value2, list):
            if len(value1) != len(value2):
                _match = False
            # check if the lists contain dictionaries as elements
            if len(value1) > 0 and isinstance(value1[0], dict_types) and isinstance(value2[0], dict_types):
                for i in range(len(value1)):
                    if not compare_dicts(value1[i], value2[i], logger):
                        _match = False
                        break
            elif sorted(value1) != sorted(value2):
                _match = False
        elif isinstance(value1, np.ndarray) and isinstance(value2, np.ndarray):
            if not np.array_equal(value1, value2):
                _match = False
        else:
            try:
                if value1 != value2:
                    _match = False
            except ValueError:
                logger.exception(type(value1))
                _match = False
        if not _match:
            logger.error(f"Value of key {key} not equal")
            match = False

    return match


def test(dataset, path: Optional[str] = None):
    logger_dataset = main_logger.getChild(f"Version_{dataset['id']}")

    def run_test(_path):
        current_labeled_h5ad_path = f"{_path}/{dataset['id']}_current_labeled.h5ad"
        raw_h5ad_path = f"{_path}/{dataset['id']}_raw.h5ad"
        new_labeled_h5ad_path = f"{_path}/{dataset['id']}_new_labeled.h5ad"

        download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path, logger_dataset)
        logger_dataset.info("Validating")
        validate(raw_h5ad_path, new_labeled_h5ad_path)

        # Compare h5ads
        logger_dataset.info("Comparing")
        new_ad = ad.read_h5ad(new_labeled_h5ad_path, backed="r")
        current_ad = ad.read_h5ad(current_labeled_h5ad_path, backed="r")
        if not current_ad.obs.equals(new_ad.obs):
            logger_dataset.error("obs not equal")
        if not current_ad.var.equals(new_ad.var):
            logger_dataset.error("var not equal")
        new_ad.uns["citation"] = current_ad.uns["citation"]
        compare_dicts(current_ad.uns, new_ad.uns, logger_dataset)

    try:
        if path:
            run_test(path)
        else:
            with tempfile.TemporaryDirectory() as path:
                run_test(path)
    except Exception as e:
        logger_dataset.exception(f"Error testing dataset: {e}")


def main():
    test_datasets_ids = dataset_ids_to_test()
    dataset_ids = "\n\t".join(test_datasets_ids)
    main_logger.info(f"Datasets Versions: {dataset_ids}")
    with logging_redirect_tqdm():
        # from tqdm import tqdm
        # for dataset in tqdm(iterate_dataset_metadata(test_datasets_ids)):
        #     test(dataset)
        from tqdm.contrib import concurrent

        concurrent.process_map(test, iterate_dataset_metadata(test_datasets_ids), max_workers=4)


if __name__ == "__main__":
    main()
