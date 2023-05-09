#!/usr/bin/env python

import gzip
import json
import os
import requests

API_URL = {
    "prod": "https://api.cellxgene.cziscience.com",
    "staging": "https://api.cellxgene.staging.single-cell.czi.technology",
    "dev": "https://api.cellxgene.dev.single-cell.czi.technology",
}


def dry_run():
    # Load processed ontologies file
    ontologies = "cellxgene_schema_cli/cellxgene_schema/ontology_files/all_ontology.json.gz"
    with gzip.open(ontologies, "r") as f:
        onto_map = json.loads(f.read().decode('utf-8'))

    # fetch all currently published datasets
    # TODO: include datasets from private collections (requires curator auth, do not upload to GHA)
    base_url = API_URL[os.getenv("env", default="dev")]
    datasets = requests.get(f"{base_url}/curation/v1/datasets").json()
    # dataset metadata fields that contain ontology terms
    ontology_types = {
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "organism",
        "self_reported_ethnicity",
        "sex",
        "tissue",
    }
    # cache terms we know are not deprecated to skip processing; init with special-case, non-ontology terms we use
    non_deprecated_term_set = {"multiethnic", "unknown", "na"}
    f = open("ontologies-curator-report.txt", "w")
    # for every dataset, check its ontology term metadata to see if any terms are deprecated. If so, report.
    for dataset in datasets:
        for ontology_type in ontology_types:
            for ontology_term in dataset[ontology_type]:
                ontology_term_id = ontology_term["ontology_term_id"]
                if ontology_term_id in non_deprecated_term_set:
                    continue
                else:
                    # all_ontologies is indexed by prefix and term ID without suffixes, so we must build the search term
                    ontology_id_parts = ontology_term_id.split(" ")[0].split(":")
                    term_prefix = ontology_id_parts[0]
                    ontology_index_id = f"{term_prefix}:{ontology_id_parts[1]}"
                    ontology = onto_map[term_prefix][ontology_index_id]

                if ontology["deprecated"]:
                    if not "replaced_by" not in ontology:
                        f.write("ALERT: Requires Manual Curator Intervention\n")
                    f.write(f"Collection ID: {dataset['collection_id']}\n")
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
                    non_deprecated_term_set.add(ontology_term_id)
    f.close()


if __name__ == "__main__":
    dry_run()
