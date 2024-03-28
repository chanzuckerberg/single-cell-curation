import hashlib
import json
import logging
import random
import tempfile
import threading
from typing import Any, Dict, List, Optional
from urllib.request import urlopen

import anndata as ad
import boto3
from cellxgene_schema.validate import validate
from tqdm import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm

logger = logging.getLogger(__name__)


def dataset_ids_to_test() -> list[str]:
    logger.info("Fetching dataset ids to test")
    # list all of the datasets
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/curation/v1/datasets").read().decode("utf-8"))
    logger.info(f"Found {len(datasets)} datasets")
    # Pick 20%
    sample_size = int(len(datasets) * 0.2)
    logger.info(f"Sampling {sample_size} datasets")
    # Only test datasets that are less than 3GB
    filtered_datasets = [
        dataset["dataset_version_id"] for dataset in datasets if dataset["assets"][0]["filesize"] < 3e9
    ]
    sample = random.sample(filtered_datasets, sample_size)

    return sample


def iterate_dataset_metadata(dataset_ids) -> List[Dict[str, Any]]:
    datasets = json.loads(urlopen("https://api.cellxgene.cziscience.com/dp/v1/datasets/index").read().decode("utf-8"))
    filtered_dataset_ids = []
    for dataset in datasets:
        if dataset["id"] in dataset_ids:
            filtered_dataset_ids.append(dataset)
    return filtered_dataset_ids


def download(uri, local_path):
    # Download raw to a temp file
    s3 = boto3.client("s3")
    logger.info(f"Downloading {uri} to {local_path}")
    bucket, key = uri.split("//", 1)[1].split("/", 1)
    s3.download_file(bucket, key, local_path)
    logger.info(f"Downloaded {uri} to {local_path}")


def download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path):
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


def hash_file(file_name: str) -> str:
    with open(file_name, "rb") as f:
        # Read the contents of the file in chunks
        chunk_size = 1024
        hasher = hashlib.sha256()
        while chunk := f.read(chunk_size):
            hasher.update(chunk)
    return hasher.hexdigest()


def test(dataset, path: Optional[str] = None):
    def run_test(_path):
        current_labeled_h5ad_path = f"{_path}/{dataset['id']}_current_labeled.h5ad"
        raw_h5ad_path = f"{_path}/{dataset['id']}_raw.h5ad"
        new_labeled_h5ad_path = f"{_path}/{dataset['id']}_new_labeled.h5ad"

        download_dataset(dataset, raw_h5ad_path, current_labeled_h5ad_path)
        logger.info(f"Validating dataset {dataset['id']}")
        validate(raw_h5ad_path, new_labeled_h5ad_path)
        logger.info(f"Comparing {dataset['id']}")

        # Compare h5ads
        new_ad = ad.read_h5ad(new_labeled_h5ad_path, backed="r")
        current_ad = ad.read_h5ad(current_labeled_h5ad_path, backed="r")
        new_ad.uns["citation"] = (current_ad.uns["citation"],)
        if not all(
            [current_ad.obs.equals(new_ad.obs), current_ad.var.equals(new_ad.var), current_ad.uns == new_ad.uns]
        ):
            logger.error(f"Datasets {dataset['id']} are not equal")

    if path:
        run_test(path)
    else:
        with tempfile.TemporaryDirectory() as path:
            run_test(path)


def main():
    logging.basicConfig(level=logging.ERROR)
    test_datasets_ids = dataset_ids_to_test()
    with logging_redirect_tqdm():
        for dataset in tqdm(iterate_dataset_metadata(test_datasets_ids)):
            logger.info(f"Testing dataset {dataset['id']}")
            test(dataset)
        #     concurrent.process_map(test, iterate_dataset_metadata(test_datasets_ids), max_workers=2)


if __name__ == "__main__":
    main()
