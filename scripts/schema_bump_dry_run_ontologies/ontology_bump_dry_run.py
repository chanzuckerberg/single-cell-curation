#!/usr/bin/env python
import json
import os
from collections import defaultdict

from scripts.common.thirdparty.discovery_api import (
    BASE_API,
    fetch_private_collections,
    fetch_private_dataset,
    fetch_public_datasets,
    get_headers,
)
from cellxgene_ontology_guide.ontology_parser import OntologyParser

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


def map_deprecated_terms(
    curator_report_entry_map: dict,  # type: ignore
    dataset: dict,  # type: ignore
    onto_parser: OntologyParser,  # type: ignore
    non_deprecated_term_cache: set,  # type: ignore
    replaced_by_map: dict,  # type: ignore
) -> None:
    """
    For a dataset, detects all deprecated ontology terms and, for each found, populates an entry in
    curator_report_entry_map with data required to report the deprecated term and any ontology-provided guidance
    to replace it. If a Deprecated Term is detected with a single, same-ontology 'replaced by' term, it is also
    added to the replaced_by_map.
    :param curator_report_entry_map: dict with information on deprecated terms detected in the data-portal corpus, and
    guidance from ontologies on how to replace them
    :param dataset: dict with dataset metadata
    :param onto_parser: OntologyParser, used to fetch ontology term metadata for a given schema_version
    :param non_deprecated_term_cache: set caching known terms that are not deprecated, used to skip re-processing
    :param replaced_by_map: dict mapping deprecated terms to a single known same-ontology replacement terms,
    updated in-place
    """
    collection_id = dataset["collection_id"]
    for ontology_type in ONTOLOGY_TYPES:
        if ontology_type in dataset:
            for ontology_term in dataset[ontology_type]:
                ontology_term_ids = ontology_term["ontology_term_id"].split(",")
                for ontology_term_id in ontology_term_ids:
                    if ontology_term_id in non_deprecated_term_cache:
                        continue
                    if onto_parser.is_term_deprecated(ontology_term_id):
                        if ontology_term_id in curator_report_entry_map[collection_id]:
                            curator_report_entry_map[collection_id][ontology_term_id]["dataset_ct"] += 1
                        else:
                            entry = dict()
                            entry["needs_alert"] = False
                            entry["dataset_ct"] = 1  # type: ignore
                            ontology_metadata = onto_parser.get_term_metadata(ontology_term_id)
                            replacement_term = onto_parser.get_term_replacement(ontology_term_id)
                            if ontology_metadata["term_tracker"]:
                                entry["term_tracker"] = ontology_metadata["term_tracker"]
                            if ontology_metadata["comments"]:
                                entry["comments"] = ontology_metadata["comments"]
                            if replacement_term:
                                entry["replaced_by"] = replacement_term
                                if ontology_metadata["replaced_by"].split(":")[0] != ontology_term_id.split(":")[0]:
                                    entry["needs_alert"] = True
                                else:
                                    if ontology_term_id not in replaced_by_map[ontology_type]:
                                        replaced_by_map[ontology_type][ontology_term_id] = replacement_term
                            else:
                                entry["needs_alert"] = True
                                if ontology_metadata["consider"]:
                                    entry["consider"] = ontology_metadata["consider"]

                            curator_report_entry_map[collection_id][ontology_term_id] = entry
                    else:
                        non_deprecated_term_cache.add(ontology_term_id)


def write_to_curator_report(output_file: str, curator_report_entry_map: dict, revision_map: dict = None) -> None:  # type: ignore
    """
    Writes curator report entries to output_file. Each entry reports a Collection ID, a Deprecated Term that has been
    detected in that collection, how many datasets are affected, and info derived from owl files that guide how to
    replace the deprecated term.
    :param output_file: filepath to write curator report to
    :param curator_report_entry_map: map containing info for curator report entries, Collection ID -> Deprecated Term
    :param revision_map: mapping of Revision ID -> Collection ID for all revisions
    """
    with open(output_file, "a") as f:
        for collection_id in curator_report_entry_map:
            for deprecated_term, entry in curator_report_entry_map[collection_id].items():
                if entry["needs_alert"]:
                    f.write("ALERT: Requires Manual Curator Intervention\n")
                f.write(f"Collection ID: {collection_id}\n")
                if revision_map and collection_id in revision_map:
                    f.write(f"Note--In A Revision of: {revision_map[collection_id]}\n")
                f.write(f"# of Affected Datasets: {entry['dataset_ct']}\n")
                f.write(f"Deprecated Term: {deprecated_term}\n")
                if "replaced_by" in entry:
                    f.write(f"Replaced By: {entry['replaced_by']}\n")
                if "consider" in entry:
                    f.write(f"Consider: {entry['consider']}\n")
                if "term_tracker" in entry:
                    f.write(f"Term Tracker: {entry['term_tracker']}\n")
                if "comments" in entry:
                    f.write(f"Comments: {entry['comments']}\n")
                f.write("\n")


def dry_run(curator_report_filepath: str, replaced_by_filepath: str) -> None:
    """
    main function which coordinates fetching information to populate a curator report and writing it to the output file,
    as well as writing a JSON file mapping any deprecated terms to known replacement terms in the same ontology.
    :param curator_report_filepath: filepath to write curator report to
    :param replaced_by_filepath: filepath to write replaced by JSON mapping to
    """
    # cache terms we know are not deprecated to skip processing; init with special-case, non-ontology terms we use
    non_deprecated_term_cache = {"unknown", "na"}
    # map deprecated terms with known, deterministic 'replaced by' terms
    replaced_by_map = {ontology_type: dict() for ontology_type in ONTOLOGY_TYPES}  # type: ignore

    base_url = BASE_API[os.getenv("corpus_env", default="dev")]
    datasets = fetch_public_datasets(base_url)  # type: ignore
    public_curator_report_entry_map = defaultdict(dict)  # type: ignore
    onto_parser = OntologyParser(schema_version="v5.0.0")
    with open(curator_report_filepath, "w") as f:
        # for every dataset, check its ontology term metadata to see if any terms are deprecated. If so, report.
        f.write("Deprecated Terms in Public Datasets:\n\n")
    for dataset in datasets:
        map_deprecated_terms(
            public_curator_report_entry_map,
            dataset,
            onto_parser,
            non_deprecated_term_cache,
            replaced_by_map,
        )
    write_to_curator_report(curator_report_filepath, public_curator_report_entry_map)

    headers = get_headers(base_url)  # type: ignore
    private_collections = fetch_private_collections(base_url, headers)  # type: ignore
    with open(curator_report_filepath, "a") as f:
        f.write("\nDeprecated Terms in Private Datasets:\n\n")
    revision_map = dict()
    private_curator_report_entry_map = defaultdict(dict)  # type: ignore
    for collection in private_collections:
        collection_id = collection["collection_id"]
        if collection.get("revision_of"):
            revision_map[collection_id] = collection["revision_of"]
        for ds in collection["datasets"]:
            dataset_id = ds["dataset_id"]
            # TODO: consider adding ontology fields to dataset preview response so a follow-up call isn't needed
            dataset_metadata = fetch_private_dataset(base_url, headers, collection_id, dataset_id)  # type: ignore
            # only process uploaded datasets
            if "processing_status" not in dataset_metadata or dataset_metadata["processing_status"] != "SUCCESS":
                continue
            map_deprecated_terms(
                private_curator_report_entry_map,
                dataset_metadata,
                onto_parser,
                non_deprecated_term_cache,
                replaced_by_map,
            )
    write_to_curator_report(curator_report_filepath, private_curator_report_entry_map, revision_map)

    if revision_map:
        with open(curator_report_filepath, "a") as f:
            f.write("\nThe Following Public Collections Will Not Be Auto-Migrated Due To Having an Open Revision:\n")
            for revision in revision_map.values():
                f.write(f"{revision}\n")

    with open(replaced_by_filepath, "w") as j:
        json.dump(replaced_by_map, j)


if __name__ == "__main__":
    dry_run("ontologies-curator-report.txt", "replaced-by.json")
