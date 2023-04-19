#!/usr/bin/env python

import requests
from owlready2 import World
from yaml import safe_load

if __name__ == "__main__":
    # Load owl.info to grab latest ontology sources
    owl_info_yml = "cellxgene_schema_cli/cellxgene_schema/ontology_files/owl_info.yml"
    with open(owl_info_yml, "r") as owl_info_handle:
        owl_info = safe_load(owl_info_handle)
    onto_map = dict()
    for onto, onto_info in owl_info.items():
        latest = onto_info["latest"]
        onto_map[onto] = World().get_ontology(onto_info["urls"][latest]).load()

    # get set of onto_term_ids for each type in our (public) corpus
    datasets = requests.get("https://api.cellxgene.dev.single-cell.czi.technology/curation/v1/datasets").json()
    ontology_names = {
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "organism",
        "self_reported_ethnicity",
        "sex",
        "tissue",
    }
    deprecated_term_map = dict()
    # seed with known special-case terms we use
    non_deprecated_term_set = {"multiethnic", "unknown", "na"}
    f = open("curator-report.txt", "w")
    for dataset in datasets:
        for key in ontology_names:
            for ontology_term in dataset[key]:
                ontology_term_id = ontology_term["ontology_term_id"]
                if ontology_term_id in non_deprecated_term_set:
                    continue
                elif ontology_term_id in deprecated_term_map:
                    ontology = deprecated_term_map[ontology_term_id]
                else:
                    ontology_id_parts = ontology_term_id.split(" ")[0].split(":")
                    prefix = ontology_id_parts[0]
                    ontology_search_id = f"{prefix}_{ontology_id_parts[1]}"
                    onto_world = onto_map[prefix]
                    ontology = onto_world.search_one(iri=f"*{ontology_search_id}")
                    if not ontology:
                        print(f"Not Found: {ontology_search_id}")
                        continue

                if ontology.deprecated:
                    if not ontology.IAO_0100001:
                        f.write("ALERT: Requires Manual Curator Intervention\n")
                    f.write(f"Collection ID: {dataset['collection_id']}\n")
                    f.write(f"Dataset ID: {dataset['dataset_id']}\n")
                    f.write(f"Deprecated Term: {ontology_term_id}\n")
                    f.write(f"Replaced By: {ontology.IAO_0100001}\n")
                    f.write(f"Consider: {ontology.consider}\n")
                    f.write(f"Comments: {ontology.comment}\n")
                    f.write("\n")
                    deprecated_term_map[ontology_term_id] = ontology
                else:
                    non_deprecated_term_set.add(ontology_term_id)
    f.close()
