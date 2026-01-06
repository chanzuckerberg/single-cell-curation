"""Annotate guidescan2 gRNA output with gene overlaps using bioframe.

This module processes guidescan2 enumerate output and annotates gRNAs with
information about which genes they overlap using genomic interval analysis.
"""

import logging
import shutil
import subprocess

import anndata as ad
import pandas as pd

from . import env
from .gencode import SupportedOrganisms
from .reference_file_manager import ReferenceFileManager

logger = logging.getLogger(__name__)

# Organisms supported for genetic_perturbations annotation per schema 7.1.0
SUPPORTED_PERTURBATION_ORGANISMS = {
    "NCBITaxon:7955",  # Danio rerio
    "NCBITaxon:9606",  # Homo sapiens
    "NCBITaxon:10090",  # Mus musculus
}


# ============================================================================
# AnnData Extraction & Validation
# ============================================================================


def _check_guidescan2_installed() -> None:
    """Check if guidescan2 is installed and available.

    :raises RuntimeError: If guidescan command is not found
    """
    logger.info("Checking for guidescan2 installation...")
    if shutil.which("guidescan") is None:
        raise RuntimeError(
            "guidescan2 not found. Please install guidescan-cli from https://github.com/pritykinlab/guidescan-cli"
        )
    logger.info("✓ guidescan2 found")


def _get_organism_from_h5ad(adata: ad.AnnData) -> str:
    """Extract organism NCBITaxon ID from h5ad.

    :param ad.AnnData adata: AnnData object with organism_ontology_term_id in uns
    :return: NCBITaxon ID (e.g., "NCBITaxon:9606")
    :rtype: str
    :raises KeyError: If organism_ontology_term_id is not in adata.uns
    :raises ValueError: If organism is not supported for genetic perturbations annotation
    """
    if "organism_ontology_term_id" not in adata.uns:
        raise KeyError("adata.uns does not contain 'organism_ontology_term_id'")

    organism_id: str = str(adata.uns["organism_ontology_term_id"])

    # Validate organism is supported for genetic perturbations
    if organism_id not in SUPPORTED_PERTURBATION_ORGANISMS:
        raise ValueError(
            f"Organism {organism_id} is not supported for genetic perturbations annotation. "
            f"Supported organisms: {', '.join(SUPPORTED_PERTURBATION_ORGANISMS)}"
        )

    # Validate organism exists in SupportedOrganisms enum
    try:
        organism = SupportedOrganisms(organism_id)
    except ValueError:
        raise ValueError(f"Invalid NCBITaxon ID: {organism_id}") from None

    logger.info(f"Detected organism: {organism.name} ({organism_id})")
    return organism_id


def _extract_perturbations_to_guidescan_csv(adata: ad.AnnData, output_csv: str) -> None:
    """Extract genetic perturbations to CSV for guidescan2 input.

    Extracts guide sequences from adata.uns['genetic_perturbations'] and writes
    them to a CSV file in guidescan2 enumerate input format.
    Guidescan enumerate requires: id, sequence, pam, chromosome, position, sense
    (chromosome, position, sense can be empty for enumeration)

    :param ad.AnnData adata: AnnData object with genetic_perturbations in uns
    :param str output_csv: Path to output CSV file
    """
    genetic_perturbations = adata.uns["genetic_perturbations"]

    # Extract guide sequences
    rows = []
    for guide_id, guide_data in genetic_perturbations.items():
        protospacer = guide_data.get("protospacer_sequence", "")
        pam_motif = guide_data.get("protospacer_adjacent_motif", "")
        pam = pam_motif[3:] if pam_motif.startswith("3' ") else pam_motif

        rows.append(
            {
                "id": guide_id,
                "sequence": protospacer,
                "pam": pam,
                "chromosome": None,
                "position": None,
                "sense": None,
            }
        )

    # Write to CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False, na_rep="")
    logger.info(f"Extracted {len(rows)} guide sequences to {output_csv}")


# ============================================================================
# Manage GuideScan2 Index Files
# ============================================================================

GUIDESCAN_INDEX_CATEGORY = "guidescan_indexes"


def _get_or_download_index(organism_id: str) -> str:
    """Get cached index or download if not present. Returns index prefix path."""
    manager = ReferenceFileManager(env.REFERENCE_CACHE_DIR, env.REFERENCE_FILES_YAML)
    key = manager.get_key_by_organism_id(GUIDESCAN_INDEX_CATEGORY, organism_id)
    if key is None:
        raise ValueError(f"GuideScan2 index not available for organism {organism_id}.")

    # pooch handles caching - returns list of extracted files
    files = manager.fetch(GUIDESCAN_INDEX_CATEGORY, key)

    # Find the index prefix from .gs file (guidescan2 expects prefix.gs, prefix.forward, prefix.reverse)
    for f in files:
        if f.endswith(".gs"):
            return f[: -len(".gs")]

    raise ValueError(f"No .gs file found in downloaded index for {organism_id}")


# ============================================================================
# GuideScan2 Integration
# ============================================================================


def _run_guidescan_enumerate(input_csv: str, index_prefix: str, output_csv: str) -> None:
    """Run guidescan enumerate command.

    :param str input_csv: Path to input CSV with guide sequences
    :param str index_prefix: Path prefix to guidescan2 index files
    :param str output_csv: Path to output CSV file
    :raises RuntimeError: If guidescan enumerate fails
    """
    cmd = [
        "guidescan",
        "enumerate",
        index_prefix,
        "-f",
        input_csv,
        "-o",
        output_csv,
        "--format",
        "csv",
        "--mismatches",
        "0",
    ]

    logger.info(f"Running guidescan enumerate: {' '.join(cmd)}")

    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        logger.info("Guidescan enumerate completed successfully")
        if result.stdout:
            logger.debug(f"Guidescan stdout: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"Guidescan enumerate failed with exit code {e.returncode}")
        logger.error(f"Stdout: {e.stdout}")
        logger.error(f"Stderr: {e.stderr}")
        raise RuntimeError(f"Guidescan enumerate failed: {e.stderr}") from e
