#!/usr/bin/env python
import logging
import os
from typing import Any, Dict, List, Optional, Tuple

import tiledb as tiledb
from jinja2 import Template

from cellxgene_schema_cli.cellxgene_schema.env import ONTOLOGY_DIR
from cellxgene_schema_cli.cellxgene_schema.ontology import SupportedOrganisms
from schema_bump_dry_run_scripts.common import (
    BASE_API,
    fetch_private_collections,
    fetch_private_dataset,
    fetch_public_datasets,
    get_headers,
)
from schema_bump_dry_run_scripts.logger import configure_logging

configure_logging()
logger = logging.getLogger()

ctx = tiledb.default_ctx({"vfs.s3.region": "us-west-2"})


def generate_report(data) -> str:
    report = """## Deprecated Terms in Public Datasets:

{% for collection in deprecated_public %}
Collection ID: {{ collection }}
Number of Affected Datasets: {{ deprecated_public[collection].num_datasets }}
Number of Deprecated Terms in Collection: {{ deprecated_public[collection].num_deprecated_genes }}
Number of Terms in Collection: {{ deprecated_public[collection].num_genes }}
Deprecated Terms: 
    {{ deprecated_public[collection].deprecated_terms  | join(', ') | wordwrap(78) | replace('\n', '\n    ')}}

{% endfor %}
## Deprecated Genes in Private Datasets:

{% for collection in open_revisions %}
Collection ID: {{ collection}}
{% if open_revisions[collection].revision_of %}
Note--In A Revision of: {{ open_revisions[collection].revision_of }}
{% endif %}
Number of Affected Datasets: {{ open_revisions[collection].num_datasets }}
Number of Deprecated Terms in Collection: {{ open_revisions[collection].num_deprecated_genes }}
Number of Terms in Collection: {{ open_revisions[collection].num_genes }}
Deprecated Terms: 
    {{ open_revisions[collection].deprecated_terms | join(', ') | wordwrap(78) | replace('\n', '\n    ')}}

{% endfor %}
## The Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:
{% for collection in non_auto_migrated %}
{{ collection }}
{% endfor %}
"""

    j2_template = Template(report, trim_blocks=True, lstrip_blocks=True)
    report = j2_template.render(data)
    return report


def get_genes(dataset: dict, stage: str) -> List[str]:
    """
    Uses tiledb to get the genes for a dataset. This method is slower, but does not add a dependency on the explorer. It
    is also free if we run computer with in the same AWS region.

    :param dataset: dataset metadata
    :param stage: prod, staging, or dev
    :return:
    """
    dataset_version_id = dataset["dataset_version_id"]
    s3_path = (
        f"s3://hosted-cellxgene-{stage}/{dataset_version_id}.cxg/var"
        if dataset.get("s3_uri")
        else dataset["s3_uri"] + "/var"
    )
    with tiledb.open(s3_path, "r") as var:
        var_df = var.df[:]
        suffix = 0
        while f"name_{suffix}" not in var_df.columns:
            suffix += 1
        index_name = "name_{suffix}"
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
            # only process uploaded datasets
            if "processing_status" not in dataset_metadata or dataset_metadata["processing_status"] != "SUCCESS":
                continue
            yield dataset_metadata, collection.get("revision_of")


def compare_genes(
    dataset: Dict[str, Any], diff_map: Dict[str, str], deprecated_datasets: Dict[str, Dict[str, str]]
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
    dataset_id = dataset["dataset_id"]
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
        if collection_id not in deprecated_datasets:
            deprecated_datasets[collection_id] = {"num_datasets": 0, "deprecated_genes": set(), "genes": set()}
        collection = deprecated_datasets[collection_id]
        collection["genes"].update(dataset_genes_to_compare)
        collection["deprecated_genes"].update(deprecated_genes_in_dataset)
        collection["num_datasets"] += 1
        logger.info(f"Dataset {dataset_id} has {len(deprecated_genes_in_dataset)} deprecated genes")
    else:
        logger.info(f"Dataset {dataset_id} has no deprecated genes")

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
    datasets = fetch_public_datasets(base_url)
    public_deprecated = {}
    for dataset in datasets:
        public_deprecated, _ = compare_genes(dataset, diff_map, public_deprecated)
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
    private_deprecated = {}
    non_auto_migrated = []
    for dataset, revision_of in fetch_private_datasets(base_url):
        private_deprecated, is_deprecated_genes_found = compare_genes(dataset, diff_map, private_deprecated)
        if revision_of and revision_of not in non_auto_migrated and is_deprecated_genes_found:
            non_auto_migrated.append(revision_of)
            private_deprecated[dataset["collection_id"]]["revision_of"] = revision_of
    return private_deprecated, non_auto_migrated


def main():
    base_url = base_url = BASE_API[os.getenv("corpus_env", default="dev")]

    report_data = {}
    diff_map = get_diff_map()

    report_data["deprecated_public"] = generate_deprecated_public(base_url, diff_map)
    report_data["open_revisions"], report_data["non_auto_migrated"] = generate_deprecated_private(base_url, diff_map)

    report = generate_report(report_data)
    with open("genes-curator-report.txt", "w") as fp:
        fp.write(report)


if __name__ == "__main__":
    main()
