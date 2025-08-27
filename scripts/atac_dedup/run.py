#!/usr/bin/env python3
"""
This script is used to run the ATAC-seq deduplication process.

Steps:
- use the api key in /Users/trentsmith/workspace/single-cell-curation/notebooks/curation_api/api-dev.key to
authenticate to the curation api
- Then use the curation api to get the list of collections that have atac-seq data. The python code is in
/Users/trentsmith/workspace/single-cell-curation/notebooks/curation_api/python/src/
- find all collection with ATAC-seq data
- track which collections are public and which are private
- create a revision for all public collections
- for each atac-seq dataset in each collection,
 - get the manifest
 - download the fragment file
 - deduplicate the fragment file
 - upload the deduplicated fragment file to S3
 - update the manifest to point to the new deduplicated fragment file
 - submit the updated manifest
- any collection that was public should be publish. All private collections should remain private.
"""
import json
import logging
import shutil
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add the curation API src to the path
curation_api_path = Path(__file__).parent.parent.parent / "notebooks/curation_api/python"
sys.path.append(str(curation_api_path))

from src.collection import create_revision, get_collections
from src.dataset import (
    get_dataset_manifest,
    upload_datafiles_from_manifest,
)
from src.utils.config import set_api_access_config

# Add the cellxgene_schema to the path
schema_path = Path(__file__).parent.parent.parent / "cellxgene_schema_cli"
sys.path.append(str(schema_path))

from cellxgene_schema import atac_seq

# Setup logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def find_atac_collections() -> Tuple[List[Dict], List[Dict]]:
    """
    Find all collections that contain ATAC-seq datasets.

    Returns:
        Tuple of (public_collections, private_collections)
    """
    logger.info("Fetching all collections...")

    # Get both public and private collections
    # INSERT_YOUR_CODE
    import concurrent.futures

    def _get_collections_with_visibility(visibility=None):
        if visibility is None:
            return get_collections()
        else:
            return get_collections(visibility=visibility)

    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_public = executor.submit(_get_collections_with_visibility, None)
        future_private = executor.submit(_get_collections_with_visibility, "PRIVATE")
        public_collections = future_public.result()
        private_collections = future_private.result()

    all_collections = public_collections + private_collections

    logger.info(f"Found {len(all_collections)} total collections")

    atac_public_collections = []
    atac_private_collections = []

    for collection in all_collections:
        has_atac = False
        for dataset in collection.get("datasets", []):
            assay_ontology_terms = [a["ontology_term_id"] for a in dataset.get("assay", [])]
            if any(atac_seq.is_atac(term) != "n" for term in assay_ontology_terms):
                has_atac = True
                break

        if has_atac:
            if collection["visibility"] == "PUBLIC":
                atac_public_collections.append(collection)
            else:
                atac_private_collections.append(collection)

    logger.info(f"Found {len(atac_public_collections)} public collections with ATAC-seq data")
    logger.info(f"Found {len(atac_private_collections)} private collections with ATAC-seq data")

    return atac_public_collections, atac_private_collections


def get_atac_datasets_from_collection(collection: Dict) -> List[Dict]:
    """Get all ATAC-seq datasets from a collection."""
    atac_datasets = []

    for dataset in collection.get("datasets", []):
        collection_id = collection["collection_id"]
        dataset_id = dataset["dataset_id"]
        manifest = get_dataset_manifest(collection_id, dataset_id)
        if manifest.get("atac_fragment"):
            atac_datasets.append(dataset)
        # assay_ontology_terms = [a["ontology_term_id"] for a in dataset.get("assay", [])]
        # if any(atac_seq.is_atac(term) != "n" for term in assay_ontology_terms):
        #     atac_datasets.append(dataset)

    return atac_datasets


ATAC_SIZE = []
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
        manifest = get_dataset_manifest(collection_id, dataset_id)
        logger.info(f"Retrieved manifest with {len(manifest)} files")

        # Look for fragment files in the manifest
        fragment_files = {k: v for k, v in manifest.items() if "fragment" in k.lower()}
        if not fragment_files:
            return None
        manifest["flags"] = {"deduplicate_fragments": True}
        upload_datafiles_from_manifest(collection_id, dataset_id, manifest)
        return True
    except Exception as e:
        logger.exception(f"Error processing dataset {dataset_id}: {e}")
        return False


def main():
    """Main execution function."""
    if shutil.which("sort") is None:
        raise RuntimeError(
            "The 'sort' command is not installed or not found in PATH. It is required to sort the fragment file."
        )
    if shutil.which("bgzip") is None:
        raise RuntimeError(
            "The 'bgzip' command is not installed or not found in PATH. It is required to compress the fragment file."
        )
    if shutil.which("pigz") is None:
        raise RuntimeError(
            "The 'pigz' command is not installed or not found in PATH. It is required to compress the fragment file."
        )
    # Configuration
    api_key_path = "./api-dev.key"

    logger.info("Starting ATAC-seq deduplication process")

    try:
        # Setup API access
        set_api_access_config(api_key_path, "staging")

        # Find collections with ATAC-seq data
        public_collections, private_collections = find_atac_collections()

        # Create revisions for public collections
        with open("./public_to_revision_map.json") as json_file:
            public_revision_mapping = json.load(json_file)
        for collection in public_collections:
            if collection["collection_id"] in public_revision_mapping:
                logger.info(f"Revision already exists for collection {collection['collection_id']}, skipping creation")
                continue
            collection_id = collection["collection_id"]
            try:
                revision_id = create_revision(collection_id)
                public_revision_mapping[collection_id] = revision_id
                logger.info(f"Created revision {revision_id} for public collection {collection_id}")
            except Exception as e:
                logger.exception(f"Failed to create revision for collection {collection_id}: {e}")
                continue
        with open("./public_to_revision_map.json", "w") as json_file:
            json.dump(public_revision_mapping, json_file)

        # Process all collections
        all_collections = public_collections + private_collections
        successful_collections = []
        failed_collections = []

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
            atac_datasets = get_atac_datasets_from_collection(collection)

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
        logger.info(f"ATAC_SIZE: {max(ATAC_SIZE), min(ATAC_SIZE), sum(ATAC_SIZE) / len(ATAC_SIZE)}")

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


if __name__ == "__main__":
    main()
