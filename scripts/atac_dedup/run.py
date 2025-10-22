#!/usr/bin/env python3
"""
# Description
This script is used to run the ATAC-seq deduplication process. It identifies collections with ATAC-seq datasets,
downloads the fragment files, deduplicates them, and re-uploads the deduplicated files back to the server.
It handles both public and private collections, creating revisions for public collections as needed. Publishing must be
done manually using the UI.

# Setup
To run this script you need to have access to the curation API and cellxgene_schema packages. You will also need an
API key with the appropriate permissions.

## Curation API:
https://github.com/chanzuckerberg/cellxgene-curation/tree/main/notebooks/curation_api/python

## Cellxgene Schema CLI:
https://github.com/chanzuckerberg/cellxgene-schema-cli

## API Key:
An API key can be generated from the Data Portal UI under your user settings.

This script relies on the reverted changes made in: https://github.com/chanzuckerberg/single-cell-data-portal/pull/7678

"""
import concurrent.futures
import datetime
import json
import logging
import sys
from pathlib import Path
from typing import Dict, List, Optional

# Add the curation API src to the path
curation_api_path = Path(__file__).parent.parent.parent / "notebooks/curation_api/python"
sys.path.append(str(curation_api_path))

from src.collection import create_revision, get_collection, get_collections
from src.dataset import get_dataset_manifest, upload_datafiles_from_manifest
from src.utils.config import set_api_access_config

# Add the cellxgene_schema to the path
schema_path = Path(__file__).parent.parent.parent / "cellxgene_schema_cli"
sys.path.append(str(schema_path))

# Setup logging
# output the logs to a file with date and time in the filename
logging.basicConfig(
    filename=f"atac_deduplication_{datetime.datetime.now().strftime('%Y%m%d_%H%M%S')}.log",
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def find_atac_collections(collection_ids: List[str]) -> List[Dict]:
    """
    Find all collections with ATAC-seq datasets, fetching metadata and manifests concurrently.

    Args:
        collection_ids: List of collection IDs

    Returns:
        List of collections with ATAC-seq datasets (with updated manifests)
    """

    collections_with_atac = []

    def fetch_manifest_for_dataset(collection_id, dataset):
        manifest = get_dataset_manifest(collection_id, dataset["dataset_id"])
        manifest.setdefault("flags", {})["deduplicate_fragments"] = True
        dataset["manifest"] = manifest
        return dataset

    def fetch_collection_and_manifests(collection_id):
        collection = get_collection(collection_id)
        atac_datasets = [
            dataset
            for dataset in collection.get("datasets", [])
            if any(asset["filetype"] == "ATAC_FRAGMENT" for asset in dataset.get("assets", []))
        ]
        if not atac_datasets:
            return None

        # Fetch manifests concurrently for ATAC_FRAGMENT datasets
        with concurrent.futures.ThreadPoolExecutor() as manifest_executor:
            futures = [
                manifest_executor.submit(fetch_manifest_for_dataset, collection_id, dataset)
                for dataset in atac_datasets
            ]
            updated_datasets = [future.result() for future in concurrent.futures.as_completed(futures)]

        collection["datasets"] = updated_datasets
        return collection

    # Fetch collections concurrently
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(fetch_collection_and_manifests, cid) for cid in collection_ids]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                logger.error(f"Error fetching collection: {e}")
                result = None
            if result:
                collections_with_atac.append(result)

    return collections_with_atac


def get_atac_datasets_from_collection(collection: Dict) -> List[Dict]:
    """Get all ATAC-seq datasets from a collection."""
    atac_datasets = []
    collection_resp = get_collection(collection["collection_id"])
    for dataset in collection_resp.get("datasets", []):
        if [asset for asset in dataset.get("assets", []) if asset["filetype"] == "ATAC_FRAGMENT"]:
            dataset["collection_id"] = collection["collection_id"]
            atac_datasets.append(dataset)

    return atac_datasets


GB = 1024 * 1024 * 1024


def process_dataset(collection_id: str, dataset: Dict) -> Optional[bool]:
    """
    Process a single ATAC-seq dataset by downloading, deduplicating, and re-uploading fragment files.

    Args:
        collection_id: Collection ID
        dataset: Dataset dictionary
        work_dir: Working directory for temporary files

    Returns:
        True if processing was successful, False otherwise
    """
    dataset_id = dataset["dataset_id"]
    logger.info(f"Processing dataset {dataset_id} in collection {collection_id}")

    try:
        # Get the current manifest
        manifest = dataset["manifest"]

        # Look for fragment files in the manifest
        upload_datafiles_from_manifest(manifest, collection_id, dataset_id)
        return True
    except Exception as e:
        logger.exception(f"Error processing dataset {dataset_id}: {e}")
        return False


def main():
    """Main execution function."""
    FIND_ATAC_COLLECTIONS = False
    CREATE_REVISIONS = False
    RUN_DEDUPLICATION = True
    # Configuration
    api_key_path = "./api-prod.key"

    logger.info("Starting ATAC-seq deduplication process")

    try:
        # Setup API access
        set_api_access_config(api_key_path)

        if FIND_ATAC_COLLECTIONS:

            try:
                with open("./public_collections.json") as json_file:
                    _public_collections = json.load(json_file)
            except FileNotFoundError:
                _public_collections = get_collections()
                with open("./public_collections.json", "w") as json_file:
                    json.dump(_public_collections, json_file)
            public_collection_ids = [col["collection_id"] for col in _public_collections]
            public_atac_collections = find_atac_collections(public_collection_ids)
            logger.info(f"Found {len(public_atac_collections)} public collections with ATAC-seq data")

            try:
                with open("./private_collections.json") as json_file:
                    _private_collections = json.load(json_file)
            except FileNotFoundError:
                _private_collections = get_collections(visibility="PRIVATE")
                with open("./private_collections.json", "w") as json_file:
                    json.dump(_private_collections, json_file)
            private_collection_ids = [col["collection_id"] for col in _private_collections]
            private_atac_collections = find_atac_collections(private_collection_ids)
            logger.info(f"Found {len(private_atac_collections)} private collections with ATAC-seq data")

            # save the collections to a json file
            with open("./atac_collections.json", "w") as json_file:
                json.dump({"public": public_atac_collections, "private": private_atac_collections}, json_file)

        with open("./atac_collections.json") as json_file:
            collections = json.load(json_file)
        public_collections = collections["public"]
        private_collections = collections["private"]

        if CREATE_REVISIONS:
            # Create revisions for public collections
            try:
                with open("./public_to_revision_map.json") as json_file:
                    public_revision_mapping = json.load(json_file)
            except FileNotFoundError:
                public_revision_mapping = {}

            for collection in public_collections:
                if collection["collection_id"] in public_revision_mapping:
                    logger.info(
                        f"Revision already exists for collection {collection['collection_id']}, skipping creation"
                    )
                    continue
                collection_id = collection["collection_id"]
                try:
                    revision_id = create_revision(collection_id)
                    public_revision_mapping[collection_id] = revision_id
                    collection["revision_id"] = revision_id
                    logger.info(f"Created revision {revision_id} for public collection {collection_id}")
                except Exception as e:
                    logger.exception(f"Failed to create revision for collection {collection_id}: {e}")
                    continue
            with open("./public_to_revision_map.json", "w") as json_file:
                json.dump(public_revision_mapping, json_file)

        if RUN_DEDUPLICATION:
            with open("./public_to_revision_map.json") as json_file:
                public_revision_mapping = json.load(json_file)
            # Process all collections
            all_collections = public_collections + private_collections
            successful_collections = []
            failed_collections = []

            # for collection in [collections["private"][3]]:
            for collection in all_collections:
                collection_id = collection["collection_id"]
                is_public = collection["visibility"] == "PUBLIC"

                # Use revision ID for public collections
                if is_public and collection_id in public_revision_mapping:
                    working_collection_id = public_revision_mapping[collection_id]
                    logger.info(f"Processing public collection {collection_id} via revision {working_collection_id}")
                else:
                    working_collection_id = collection_id
                    logger.info(f"Processing private collection {collection_id}")

                # Get ATAC datasets from this collection
                atac_datasets = [dataset for dataset in collection["datasets"] if dataset.get("manifest")]

                if not atac_datasets:
                    continue

                logger.info(f"Processing {len(atac_datasets)} ATAC datasets in collection {collection_id}")

                collection_success = True
                for dataset in atac_datasets:
                    success = process_dataset(working_collection_id, dataset)
                    if success is False:
                        collection_success = False
                    elif success is None:
                        logger.info(f"No fragment files found in dataset {dataset['dataset_id']}, skipping")

                if collection_success:
                    successful_collections.append(collection_id)
                    logger.info(f"Successfully processed collection {collection_id}")
                else:
                    failed_collections.append(collection_id)
                    logger.error(f"Failed to process collection {collection_id}")

            # Summary
            logger.info("ATAC-seq deduplication process completed")
            logger.info(f"Successfully processed collections: {len(successful_collections)}")
            logger.info(f"Failed collections: {len(failed_collections)}")

            # Print a link to the collections to be published
            logger.info("Links to processed public collections:")
            for collection_id in successful_collections:
                if collection_id in public_revision_mapping:
                    revision_id = public_revision_mapping[collection_id]
                    logger.info(f"https://cellxgene.cziscience.com/collections/{revision_id}")

            logger.info("Links to processed private collections:")
            for collection_id in successful_collections:
                if collection_id not in public_revision_mapping:
                    logger.info(f"https://cellxgene.cziscience.com/collections/{collection_id}")

            if failed_collections:
                logger.error(f"Failed collection IDs: {failed_collections}")
                sys.exit(1)

    except Exception as e:
        logger.exception(f"Script failed with error: {e}")
        sys.exit(1)


def get_smallest_fragment_file():
    # Configuration
    api_key_path = "./api-staging.key"

    logger.info("Starting ATAC-seq deduplication process")

    # Setup API access
    set_api_access_config(api_key_path, "staging")

    # Find collections with ATAC-seq data
    with open("./atac_collections.json") as json_file:
        collections = json.load(json_file)
    public_collections = collections["public"]
    private_collections = collections["private"]
    # public_collections, private_collections = find_atac_collections()

    # Process all collections
    all_collections = public_collections + private_collections
    atac_datasets = []
    for collection in all_collections:
        _atac_datasets = get_atac_datasets_from_collection(collection)

        if not _atac_datasets:
            continue
        atac_datasets.extend(_atac_datasets)

    # size of the fragment by manifest
    fragment_size_by_manifest = {}
    for dataset in atac_datasets:
        fragment_size = [
            asset["filesize"] // GB for asset in dataset["assets"] if asset["filetype"] == "ATAC_FRAGMENT"
        ][0]
        manifest = get_dataset_manifest(dataset["collection_id"], dataset["dataset_id"])
        dataset["manifest"] = manifest
        fragment_size_by_manifest[fragment_size] = dataset
    with open("fragment_size_by_manifest.json", "w") as f:
        json.dump(fragment_size_by_manifest, f, sort_keys=True, indent=4)
    print(fragment_size_by_manifest[min(fragment_size_by_manifest.keys())])


if __name__ == "__main__":
    # get_smallest_fragment_file()
    main()
