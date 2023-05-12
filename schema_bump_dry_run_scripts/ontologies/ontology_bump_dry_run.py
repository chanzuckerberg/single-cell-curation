#!/usr/bin/env python
import gzip
import json
import os

import requests

BASE_API = {
    "prod": "api.cellxgene.cziscience.com",
    "staging": "api.cellxgene.staging.single-cell.czi.technology",
    "dev": "api.cellxgene.dev.single-cell.czi.technology",
}

# dataset metadata fields that contain ontology terms
ONTOLOGY_TYPES = {
    "assay",
    "cell_type",
    "development_stage",
    "disease",
    "organism",
    "self_reported_ethnicity",
    "sex",
    "tissue",
}


def load_ontology_map():
    # Load processed ontologies file
    ontologies = "cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json.gz"
    with gzip.open(ontologies, "rt") as f:
        onto_map = json.loads(f.read())
    return onto_map


def fetch_public_datasets(base_url):
    return requests.get(f"https://{base_url}/curation/v1/datasets").json()


def fetch_private_collections(base_url, headers):
    return requests.get(f"https://{base_url}/curation/v1/collections?visibility=PRIVATE", headers=headers).json()


def fetch_private_dataset(base_url, headers, collection_id, dataset_id):
    return requests.get(
        f"https://{base_url}/curation/v1/collections/{collection_id}/datasets/{dataset_id}", headers=headers
    ).json()


def get_headers(base_url):
    auth0_secrets = json.loads(os.getenv("AUTH0_SECRETS"))
    response = requests.post(
        f"https://{base_url}/curation/v1/auth/token", headers={"x-api-key": f"{auth0_secrets['super_curator_api_key']}"}
    )
    access_token = response.json()["access_token"]
    return {"Authorization": f"Bearer {access_token}"}


def report_deprecated_terms(
    output_file, dataset, collection_id, onto_map, non_deprecated_term_cache, revision_of=False
):
    with open(output_file, "a") as f:
        for ontology_type in ONTOLOGY_TYPES:
            if ontology_type in dataset:
                for ontology_term in dataset[ontology_type]:
                    ontology_term_id = ontology_term["ontology_term_id"]
                    if ontology_term_id in non_deprecated_term_cache:
                        continue
                    else:
                        # all_ontologies is indexed by ontology prefix and term ID without suffixes,
                        # so we must parse + build the search term
                        ontology_id_parts = ontology_term_id.split(" ")[0].split(":")
                        term_prefix = ontology_id_parts[0]
                        ontology_index_id = f"{term_prefix}:{ontology_id_parts[1]}"
                        ontology = onto_map[term_prefix][ontology_index_id]

                    if ontology["deprecated"]:
                        replaced_in_diff_ontology = False
                        if "replaced_by" in ontology:
                            replacement_term_ontology = ontology["replaced_by"].split("_")[0]
                            if replacement_term_ontology != term_prefix:
                                replaced_in_diff_ontology = True
                        if "replaced_by" not in ontology or replaced_in_diff_ontology:
                            f.write("ALERT: Requires Manual Curator Intervention\n")
                        if revision_of:
                            f.write(f"Note--In A Revision of: {revision_of}\n")
                        f.write(f"Collection ID: {collection_id}\n")
                        f.write(f"Dataset ID: {dataset['dataset_id']}\n")
                        f.write(f"Deprecated Term: {ontology_term_id}\n")
                        if "replaced_by" in ontology:
                            f.write(f"Replaced By: {ontology['replaced_by']}\n")
                        if "consider" in ontology:
                            f.write(f"Consider: {ontology['consider']}\n")
                        if "comments" in ontology:
                            f.write(f"Comments: {ontology['comments']}\n")
                        f.write("\n")
                    else:
                        non_deprecated_term_cache.add(ontology_term_id)


def dry_run(output_file):
    onto_map = load_ontology_map()

    # cache terms we know are not deprecated to skip processing; init with special-case, non-ontology terms we use
    non_deprecated_term_cache = {"multiethnic", "unknown", "na"}

    base_url = BASE_API[os.getenv("corpus_env", default="dev")]
    datasets = fetch_public_datasets(base_url)
    with open(output_file, "w") as f:
        # for every dataset, check its ontology term metadata to see if any terms are deprecated. If so, report.
        f.write("Deprecated Terms in Public Datasets:\n\n")
    for dataset in datasets:
        report_deprecated_terms(output_file, dataset, dataset["collection_id"], onto_map, non_deprecated_term_cache)

    headers = get_headers(base_url)
    private_collections = fetch_private_collections(base_url, headers)
    open_revision_set = set()
    with open(output_file, "a") as f:
        f.write("\nDeprecated Terms in Private Datasets:\n\n")
    for collection in private_collections:
        revision_of = collection["revision_of"]
        if revision_of:
            open_revision_set.add(revision_of)
        collection_id = collection["collection_id"]
        for ds in collection["datasets"]:
            dataset_id = ds["dataset_id"]
            # TODO: consider adding ontology fields to dataset preview response so a follow-up call isn't needed
            dataset_metadata = fetch_private_dataset(base_url, headers, collection_id, dataset_id)
            # only process uploaded datasets
            if "processing_status" not in dataset_metadata or dataset_metadata["processing_status"] != "SUCCESS":
                continue
            report_deprecated_terms(
                output_file, dataset_metadata, collection_id, onto_map, non_deprecated_term_cache, revision_of
            )
    if open_revision_set:
        with open(output_file, "a") as f:
            f.write("\nThe Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:\n")
            for revision in open_revision_set:
                f.write(f"{revision}\n")


if __name__ == "__main__":
    dry_run("ontologies-curator-report.txt")
