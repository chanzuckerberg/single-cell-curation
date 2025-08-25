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
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Tuple

import requests

# Add the curation API src to the path
curation_api_path = Path(__file__).parent.parent.parent / "notebooks/curation_api/python"
sys.path.append(str(curation_api_path))

from src.collection import create_revision, get_collections
from src.dataset import (
    download_assets_from_manifest,
    get_dataset_manifest,
    upload_datafiles_from_manifest,
    upload_local_datafile,
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
    public_collections = get_collections()
    private_collections = get_collections(visibility="PRIVATE")
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
        assay_ontology_terms = [a["ontology_term_id"] for a in dataset.get("assay", [])]
        if any(atac_seq.is_atac(term) != "n" for term in assay_ontology_terms):
            atac_datasets.append(dataset)

    return atac_datasets


def process_fragment_file(fragment_file_path: str, output_dir: str) -> str:
    """
    Process (deduplicate and sort) a fragment file.

    Args:
        fragment_file_path: Path to the downloaded fragment file
        output_dir: Directory to write the processed file

    Returns:
        Path to the processed fragment file
    """

    logger.info(f"Processing fragment file: {fragment_file_path}")

    output_filename = Path(fragment_file_path).stem.replace(".tsv", "") + "_dedup.tsv.bgz"
    bgzip_output_file = Path(output_dir) / output_filename

    # Use the existing atac_seq processing functions to deduplicate and sort
    try:
        # The prepare_fragment function handles sorting and bgzip compression
        # It expects a gzipped input file and outputs a bgzipped file
        num_cores = os.cpu_count()

        bgzip_proc = subprocess.Popen(
            f"""gzip -dc "{fragment_file_path}" \
| LC_ALL=C sort -t $'\t' -k1,1 -k2,2n -k3,3n -k4,4 \
    -S 60% --parallel="{num_cores}" \
    --compress-program="pigz -p {num_cores}" \
| uniq | bgzip -@ "{num_cores}" -c > "{bgzip_output_file}" 
            """,
            shell=True,
        )
        bgzip_proc.wait()
        # with open(bgzip_output_file, "wb") as out_f:
        #     gzip_proc = subprocess.Popen(["gzip", "-dc", fragment_file_path], stdout=subprocess.PIPE)
        #     sort_proc = subprocess.Popen(
        #         [
        #             "sort",
        #             f"--parallel={num_cores}",
        #             "-t",
        #             "\t",
        #             "-k1,1",
        #             "-k2,2n",
        #             "-k3,3n",
        #             "-k4,4",
        #             "-S 40%",
        #             '--compress-program="pigz -p {num_cores}"',
        #         ],
        #         stdin=gzip_proc.stdout,
        #         stdout=subprocess.PIPE,
        #         env={**os.environ, "LC_ALL": "C"},
        #     )
        #     gzip_proc.stdout.close()
        #     # run it through unique to remove any duplicate rows that may have been introduced
        #     unique_proc = subprocess.Popen(
        #         ["uniq"],
        #         stdin=sort_proc.stdout,
        #         stdout=subprocess.PIPE,
        #     )
        #     sort_proc.stdout.close()
        #     bgzip_proc = subprocess.Popen(["bgzip", f"--threads={num_cores}", "-c"], stdin=unique_proc.stdout,
        #     stdout=out_f)
        #     unique_proc.stdout.close()
        #     bgzip_proc.wait()
        if bgzip_proc.returncode != 0:
            raise RuntimeError(f"bgzip compression failed with error code {bgzip_proc.returncode}")
        logger.info(f"bgzip compression completed successfully for {bgzip_output_file}")
        return str(bgzip_output_file)
    except Exception as e:
        logger.exception(f"Error processing fragment file {fragment_file_path}: {e}")
        raise


def process_dataset(collection_id: str, dataset: Dict, work_dir: str) -> bool:
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
            logger.warning(f"No fragment files found in dataset {dataset_id}")
            return True  # Not an error, just no fragment files to process

        logger.info(f"Found {len(fragment_files)} fragment files to process")

        # Create dataset-specific work directory
        dataset_work_dir = Path(work_dir) / dataset_id
        dataset_work_dir.mkdir(exist_ok=True)

        # Download fragment files
        original_cwd = os.getcwd()
        os.chdir(dataset_work_dir)

        try:
            # TODO: remove
            content_length = requests.head(fragment_files["atac_fragment"]).headers.get("Content-Length")
            if int(content_length) > 3 * 1024 * 1024 * 1024:  # 1 GB
                logger.warning(f"Skipping large file (>1GB): {fragment_files['atac_fragment']}")
                return True

            download_assets_from_manifest(fragment_files)

            # Process each fragment file
            updated_manifest = manifest.copy()

            for file_key, file_url in fragment_files.items():
                original_filename = file_url.split("/")[-1]
                original_file_path = dataset_work_dir / original_filename

                if not original_file_path.exists():
                    logger.warning(f"Downloaded file not found: {original_file_path}")
                    continue

                # Process the fragment file (deduplicate and sort)
                processed_file_path = process_fragment_file(str(original_file_path), str(dataset_work_dir))

                # Upload the processed file
                s3_uri = upload_local_datafile(processed_file_path, collection_id, dataset_id)

                # Update the manifest with the new file location
                updated_manifest[file_key] = s3_uri
                logger.info(f"Updated manifest entry {file_key}: {s3_uri}")

            # Upload the updated manifest
            upload_datafiles_from_manifest(updated_manifest, collection_id, dataset_id)
            logger.info(f"Successfully updated manifest for dataset {dataset_id}")

            return True
        except Exception as e:
            logger.exception(f"Error downloading or processing files for dataset {dataset_id}: {e}")
            return False
        finally:
            os.chdir(original_cwd)

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
    api_key_path = "/Users/trentsmith/workspace/single-cell-curation/notebooks/curation_api/api-dev.key"

    logger.info("Starting ATAC-seq deduplication process")

    try:
        # Setup API access
        set_api_access_config(api_key_path, "staging")

        # Find collections with ATAC-seq data
        public_collections, private_collections = find_atac_collections()

        # Create revisions for public collections
        with open(
            "/Users/trentsmith/workspace/single-cell-curation/scripts/atac_dedup/public_to_revision_map.json"
        ) as json_file:
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
                logger.info(f"No ATAC datasets found in collection {collection_id}")
                continue

            logger.info(f"Processing {len(atac_datasets)} ATAC datasets in collection {collection_id}")

            collection_success = True
            for dataset in atac_datasets:
                with tempfile.TemporaryDirectory(prefix="atac_dedup_") as work_dir:
                    success = process_dataset(working_collection_id, dataset, work_dir)
                    if not success:
                        collection_success = False

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


if __name__ == "__main__":
    main()
