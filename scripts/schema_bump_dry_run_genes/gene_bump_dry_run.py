#!/usr/bin/env python
import logging
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import tiledb as tiledb
from jinja2 import Template

from cellxgene_schema_cli.cellxgene_schema.env import ONTOLOGY_DIR
from cellxgene_schema_cli.cellxgene_schema.ontology import SupportedOrganisms
from scripts.common.thirdparty import (
    BASE_API,
    fetch_private_collections,
    fetch_private_dataset,
    fetch_public_datasets,
    get_headers,
)
from scripts.logger import configure_logging

configure_logging()
logger = logging.getLogger()

ctx = tiledb.default_ctx({"vfs.s3.region": "us-west-2"})


class RunReporter:
    def __init__(self):
        self.public_datasets_processed = 0
        self.public_deprecated_datasets = 0
        self.public_errored_datasets = 0
        self.private_datasets_processed = 0
        self.private_deprecated_datasets = 0
        self.private_errored_datasets = 0

    def log_report(self):
        logger.info("Run report:")
        logger.info(f"  Public datasets processed: {self.public_datasets_processed}")
        logger.info(f"  Public deprecated datasets: {self.public_deprecated_datasets}")
        logger.info(f"  Public errored datasets: {self.public_errored_datasets}")
        logger.info(f"  Private datasets processed: {self.private_datasets_processed}")
        logger.info(f"  Private deprecated datasets: {self.private_deprecated_datasets}")
        logger.info(f"  Private errored datasets: {self.private_errored_datasets}")


run_reporter = RunReporter()


def generate_report(data) -> str:
    file_path = os.path.dirname(os.path.realpath(__file__))
    with open(os.path.join(file_path, "report_template.jinja"), "r") as fp:
        report = fp.read()

    j2_template = Template(report, trim_blocks=True, lstrip_blocks=True)
    report = j2_template.render(data)
    return report


def get_genes(dataset: dict) -> List[str]:
    """
    Uses tiledb to get the genes for a dataset. This method is slower, but does not add a dependency on the explorer. It
    is also free if we run computer with in the same AWS region.

    :param dataset: dataset metadata
    :return:
    """

    # TODO: remove different paths once we decide if we are using data portal or curation API.
    if dataset.get("s3_uri"):
        s3_path = dataset.get("s3_uri") + "/var"
    elif dataset.get("dataset_assets"):
        s3_path = [asset["s3_uri"] for asset in dataset["dataset_assets"] if asset["filetype"] == "CXG"][0] + "/var"
    else:
        dataset_version_id = dataset["dataset_version_id"]
        stage = os.getenv("corpus_env", default="dev")
        s3_path = f"s3://hosted-cellxgene-{stage}/{dataset_version_id}.cxg/var"

    logger.debug(f"Fetching genes from {s3_path}")
    with tiledb.open(s3_path, "r") as var:
        var_df = var.df[:]
        prefix = "name_"
        suffix = 0
        filtered_columns = {i for i in var_df.columns if i.startswith(prefix)}
        while filtered_columns and suffix < len(var_df.columns) and f"{prefix}{suffix}" not in var_df.columns:
            suffix += 1
        if suffix == len(var_df.columns) or not filtered_columns:
            raise KeyError(f"No columns with matching prefix:'{prefix}' found in var_df")
        index_name = f"{prefix}{suffix}"
        stored_genes = var_df[index_name].to_list()
    return stored_genes


def get_diff_map() -> Dict[str, List[str]]:
    # list all of the files ending with diff.txt in the cellxgene_schema/ontology_files directory
    # for each file, open it and read the contents into a dictionary
    diff_map = {}
    suffix = "_diff.txt"
    files = os.listdir(ONTOLOGY_DIR)
    for file in files:
        if file.endswith(suffix):
            with open(f"{ONTOLOGY_DIR}/{file}") as fp:
                organism = getattr(SupportedOrganisms, file.removesuffix(suffix).upper()).value
                diff_map[organism] = fp.read().splitlines()
    logger.info("organisms with deprecated genes: %s", diff_map.keys())
    return diff_map


def fetch_private_datasets(base_url) -> Tuple[List[dict], Optional[str]]:
    """
    Fetches all private collections and parses the datasets from the response.
    :param base_url:
    :return:
    """
    auth_headers = get_headers(base_url)
    for collection in fetch_private_collections(base_url, auth_headers):
        collection_id = collection["collection_id"]
        for ds in collection["datasets"]:
            dataset_metadata = fetch_private_dataset(base_url, auth_headers, collection_id, ds["dataset_id"])
            dataset_metadata["collection_id"] = collection_id
            # only process uploaded datasets
            if "processing_status" not in dataset_metadata or dataset_metadata["processing_status"] != "SUCCESS":
                continue
            yield dataset_metadata, collection.get("revision_of")


def compare_genes(
    dataset: Dict[str, Any], diff_map: Dict[str, str], deprecated_datasets: defaultdict(list)
) -> Tuple[Dict, bool]:
    """
    Compare genes in a dataset with the provided diff map and update the deprecated_datasets dictionary.

    :param dataset: The dataset to compare genes for.
    :type dataset: dict
    :param diff_map: The map of organisms and their deprecated genes.
    :type diff_map: dict
    :param deprecated_datasets: The dictionary to store deprecated datasets.
    :type deprecated_datasets: dict

    :return: A tuple containing the updated deprecated_datasets dictionary and a flag indicating if any deprecated
    genes were found.
    :rtype: tuple
    """
    dataset_id = dataset.get("dataset_id") or dataset.get("id")
    collection_id = dataset["collection_id"]
    dataset_genes_to_compare = set(get_genes(dataset))
    organisms = [organism_info["ontology_term_id"] for organism_info in dataset["organism"]]
    is_deprecated_genes_found = False
    deprecated_genes_in_dataset = set()
    for organism in organisms:
        deprecated_genes_source = diff_map.get(organism, set())
        intersection_genes = dataset_genes_to_compare.intersection(deprecated_genes_source)
        if intersection_genes:
            deprecated_genes_in_dataset.update(intersection_genes)
            is_deprecated_genes_found = True

    if is_deprecated_genes_found:
        dataset_group_key = (*sorted(deprecated_genes_in_dataset), len(dataset_genes_to_compare))
        # dataset_group_key is a unique key used to identify a group of datasets that have the same set of deprecated
        # genes and the same number of genes. This is used to group datasets together in the report.
        if collection_id not in deprecated_datasets:
            # add a collection_id to the dictionary if it does not exist
            deprecated_datasets[collection_id] = {"dataset_groups": {}}
        collection = deprecated_datasets[collection_id]
        if dataset_group_key not in collection["dataset_groups"]:
            # add a dataset_group, this is used when generating the report.
            collection["dataset_groups"][dataset_group_key] = {
                "datasets": [],
                "num_datasets": 0,
                "deprecated_genes": deprecated_genes_in_dataset,
                "num_genes": len(dataset_genes_to_compare),
            }
        group = collection["dataset_groups"][dataset_group_key]
        group["num_datasets"] += 1
        group["datasets"].append(dataset_id)
        logger.debug(f"Dataset {dataset_id} has {len(deprecated_genes_in_dataset)} deprecated genes")
    else:
        logger.debug(f"Dataset {dataset_id} has no deprecated genes")

    return deprecated_datasets, is_deprecated_genes_found


def generate_deprecated_public(base_url: str, diff_map: Dict) -> Dict:
    """
    Generate a dictionary of deprecated datasets from public datasets.

    :param base_url: The base URL for fetching public datasets.
    :type base_url: str
    :param diff_map: The map of organisms and their deprecated genes.
    :type diff_map: dict

    :return: A dictionary of public collections with datasets containing deprecated genes.
    :rtype: dict
    """
    public_deprecated = {}
    for dataset in fetch_public_datasets(base_url):
        run_reporter.public_datasets_processed += 1
        try:
            public_deprecated, is_deprecated_genes_found = compare_genes(dataset, diff_map, public_deprecated)
            if is_deprecated_genes_found:
                run_reporter.public_deprecated_datasets += 1
        except Exception as e:
            run_reporter.public_errored_datasets += 1
            logger.error(
                f"Error processing public dataset {dataset['dataset_id']}, in collection {dataset['collection_id']}: "
                f"{e}"
            )
    for collection in public_deprecated.values():
        # convert dataset_groups from a dictionary to a list of its values.
        collection["dataset_groups"] = list(collection["dataset_groups"].values())
    return public_deprecated


def generate_deprecated_private(base_url: str, diff_map: Dict) -> Tuple[Dict, List]:
    """
    Generate a dictionary of deprecated private collections its datasets and a list of non-auto-migrated public
    collections.

    :param base_url: The base URL for fetching private datasets.
    :type base_url: str
    :param diff_map: The map of organisms and their deprecated genes.
    :type diff_map: dict
    :return: A tuple containing the dictionary of collections with datasets containing deprecated genes and the list
        of non-auto-migrated public collections.
    :rtype: tuple
    """
    private_deprecated = dict()
    non_auto_migrated = set()
    for dataset, revision_of in fetch_private_datasets(base_url):
        run_reporter.private_datasets_processed += 1
        try:
            private_deprecated, is_deprecated_genes_found = compare_genes(dataset, diff_map, private_deprecated)
            if is_deprecated_genes_found:
                run_reporter.private_deprecated_datasets += 1
                if revision_of and revision_of not in non_auto_migrated:
                    non_auto_migrated.add(revision_of)
                    private_deprecated[dataset["collection_id"]]["revision_of"] = revision_of
        except Exception as e:
            run_reporter.private_errored_datasets += 1
            logger.error(f"Error processing private dataset: {e}")
    for collection in private_deprecated.values():
        collection["dataset_groups"] = list(collection["dataset_groups"].values())
    return private_deprecated, non_auto_migrated


def main():
    base_url = BASE_API[os.getenv("corpus_env", default="dev")]
    logger.info(f"Using base URL: {base_url}")
    report_data = {}
    diff_map = get_diff_map()

    try:
        report_data["deprecated_public"] = generate_deprecated_public(base_url, diff_map)
    except Exception:
        logger.exception("Error generating deprecated public datasets report")
        report_data["deprecated_public"] = {}

    try:
        report_data["open_revisions"], report_data["non_auto_migrated"] = generate_deprecated_private(
            base_url, diff_map
        )
    except Exception:
        logger.exception("Error generating deprecated private datasets report")
        report_data["open_revisions"] = {}
        report_data["non_auto_migrated"] = []

    report = generate_report(report_data)
    with open("genes-curator-report.txt", "w") as fp:
        fp.write(report)
    logger.info("Curator Report generated")
    run_reporter.log_report()


if __name__ == "__main__":
    main()
