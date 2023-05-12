#!/usr/bin/env python
import json
import logging
import os

import requests
import tiledb as tiledb

from schema_bump_dry_run_scripts.genes.logger import logit, print_tracking

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
}

ctx = tiledb.default_ctx({"vfs.s3.region": "us-west-2"})


stage = "dev"
logger = logging.getLogger()


@logit
def get_datasets():
    # fetch all currently published datasets
    # TODO: include datasets from private collections (requires curator auth, do not upload to GHA)
    base_url = API_URL[os.getenv("corpus_env", default=stage)]
    datasets = requests.get(f"{base_url}/curation/v1/datasets").json()
    return datasets


@logit
def get_genes_with_explorer(dataset_version_id):
    schema = requests.get(
        f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/s3_uri/s3%253A%252F%252Fhosted-cellxgene-dev%252F{dataset_version_id}.cxg/api/v0.3/schema"
    ).json()

    feature_name = requests.get(
        f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/s3_uri/s3%253A%252F%252Fhosted"
        f"-cellxgene-dev%252F{dataset_version_id}.cxg/api/v0.3/annotations/var?annotation-name=feature_name"
    ).content
    feature_name = json.loads(feature_name[feature_name.rfind(b"[") : feature_name.rfind(b"]") + 1].decode("utf-8"))
    feature_references = []
    for col in schema["schema"]["annotations"]["var"]["columns"]:
        if col["name"] == "feature_reference":
            feature_references.extend(col["categories"])
    return {"feature_references": feature_references, "feature_name": feature_name}


@logit
def get_genes_with_tiledb(dataset_version_id):
    # s3_path = "https://hosted-cellxgene-dev.s3.amazonaws.com/{dataset_version_id}.cxg/var"
    s3_path = f"s3://hosted-cellxgene-{stage}/{dataset_version_id}.cxg/var"
    with tiledb.open(s3_path, "r") as var:
        var_df = var.df[:]
        stored_genes = set(var_df["feature_name"].to_numpy(dtype=str))
        species = json.loads(var.meta["cxg_schema"])["feature_reference"]["categories"]
    return {"feature_references": species, "feature_name": stored_genes}


def main():
    datasets = get_datasets()
    try:
        for i, dataset in enumerate(datasets):
            if i == 40:
                break
            dataset_version_id = dataset["dataset_version_id"]
            print("processed", dataset_version_id)
            try:
                get_genes_with_explorer(dataset_version_id)
                get_genes_with_tiledb(dataset_version_id)
            except Exception:
                logger.exception(f"failed {dataset_version_id}")
    finally:
        print_tracking()


if __name__ == "__main__":
    main()
