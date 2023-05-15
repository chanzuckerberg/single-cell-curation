#!/usr/bin/env python
import json
import logging
import os

import requests
import tiledb as tiledb

from schema_bump_dry_run_scripts.genes.logger import configure_logging, log_tracking, logit

configure_logging()
stage = "prod"
logger = logging.getLogger()

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
}

EXPLORER = {
    "prod": "https://api.cellxgene.cziscience.com/cellxgene/s3_uri/s3%253A%252F%252Fhosted-cellxgene-{"
    "stage}%252F{identifier}.cxg/api/v0.3/annotations/var?annotation-name=feature_name",
    "dev": "https://api.cellxgene.{stage}.single-cell.czi.technology/cellxgene/s3_uri/s3%253A%252F%252Fhosted"
    "-cellxgene-{stage}%252F{identifier}.cxg/api/v0.3/annotations/var?annotation-name=feature_name",
    "staging": "https://api.cellxgene.{stage}.single-cell.czi.technology/cellxgene/s3_uri/s3%253A%252F%252Fhosted"
    "-cellxgene-{stage}%252F{identifier}.cxg/api/v0.3/annotations/var?annotation-name=feature_name",
}
ctx = tiledb.default_ctx({"vfs.s3.region": "us-west-2"})


@logit
def get_datasets_curation():
    # fetch all currently published datasets
    # TODO: include datasets from private collections (requires curator auth, do not upload to GHA)
    base_url = API_URL[os.getenv("corpus_env", default=stage)]
    datasets = requests.get(f"{base_url}/curation/v1/datasets").json()
    return datasets


@logit
def get_datasets_dp():
    # fetch all currently published datasets
    # TODO: include datasets from private collections (requires curator auth, do not upload to GHA)
    base_url = API_URL[os.getenv("corpus_env", default=stage)]
    datasets = requests.get(f"{base_url}/dp/v1/datasets/index").json()
    for dataset in datasets:
        dataset["s3_uri"] = [i["s3_uri"] for i in dataset["dataset_assets"] if i["filetype"] == "CXG"][0]
    return datasets


@logit
def get_genes_with_explorer(identifier) -> list[str]:
    """
    Uses the explorer API to get the genes for a dataset. This depends on the explorer API ramaining constant. This
    method is faster, but adds a dependency on the explorer API, and costs more to extra data. The cost may be
    negiligble, more testing is needed to confirm.
    :param identifier:
    :return:
    """
    # This used explorer to get the speicies, but this information can be found in the dataset metadata.
    # schema = requests.get(
    #     f"https://api.cellxgene.dev.single-cell.czi.technology/cellxgene/s3_uri/s3%253A%252F%252Fhosted-cellxgene-dev%252F{identifier}.cxg/api/v0.3/schema"
    # ).json()
    # try:
    #     feature_references = []
    #     for col in schema["schema"]["annotations"]["var"]["columns"]:
    #         if col["name"] == "feature_reference":
    #             feature_references.extend(col["categories"])
    # except KeyError as ex:
    #     print(schema["message"])
    #     raise ex

    feature_name = requests.get(EXPLORER[stage].format(stage=stage, identifier=identifier)).content
    try:
        feature_name = feature_name[feature_name.rfind(b"[") : feature_name.rfind(b"]") + 1].decode("utf-8")
        feature_name = json.loads(feature_name)
    except json.decoder.JSONDecodeError as ex:
        raise ex
    return feature_name


@logit
def get_genes_with_tiledb(dataset_version_id=None, s3_uri=None) -> list[str]:
    """
    Uses tiledb to get the genes for a dataset. This method is slower, but does not add a dependency on the explorer. It
    is also free if we run computer with in the same AWS region.

    10x slower than the explorer method.

    :param dataset_version_id:
    :return:
    """
    s3_path = s3_uri + "/var" if s3_uri else f"s3://hosted-cellxgene-{stage}/{dataset_version_id}.cxg/var"
    with tiledb.open(s3_path, "r") as var:
        var_df = var.df[:]
        stored_genes = set(var_df["feature_name"].to_numpy(dtype=str))
        # species = json.loads(var.meta["cxg_schema"])["feature_reference"]["categories"] # get the speicies from tiledb
    return stored_genes


def main():
    datasets = get_datasets_dp()
    failed = []
    try:
        for i, dataset in enumerate(datasets):
            if i == 100:
                break
            # Thess calls can be run in parallel becuase they rely on network calls.
            identifiers = [dataset.get("dataset_version_id"), dataset.get("dataset_id")]
            [i["ontology_term_id"] for i in dataset["organism"]]
            try:
                # get_genes_with_explorer(identifiers[0])
                get_genes_with_tiledb(identifiers, s3_uri=dataset["s3_uri"])
            except Exception:
                failed_string = (
                    f"failed "
                    f"dataset_version_id={dataset.get('dataset_version_id')}, "
                    f"dataset_id={dataset.get('dataset_id')}, "
                    f"collection_id={dataset.get('collection_id')}, "
                    f"explorer_url={dataset.get('explorer_url')}, "
                    f"title={dataset.get('title')}"
                )
                failed.append(failed_string)
                logger.exception(failed_string)
            else:
                logger.info(f"success {identifiers[0]}")
    finally:
        if failed:
            logging.info("failed Cases")
            for fail in failed:
                logging.info(f"\t{fail}")
        log_tracking()


if __name__ == "__main__":
    main()
