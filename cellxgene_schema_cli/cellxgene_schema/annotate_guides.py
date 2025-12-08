"""Annotate guidescan2 gRNA output with gene overlaps using bioframe.

This module processes guidescan2 enumerate output and annotates gRNAs with
information about which genes they overlap using genomic interval analysis.
"""

import logging
import os
import shutil
import subprocess
import tarfile
import tempfile
import urllib.request
from typing import Dict, List, Optional, Tuple

import anndata as ad
import bioframe as bf
import pandas as pd

from . import env
from .gencode import SupportedOrganisms

logger = logging.getLogger(__name__)

# Organisms supported for genetic_perturbations annotation per schema 7.1.0
SUPPORTED_PERTURBATION_ORGANISMS = {
    "NCBITaxon:7955",  # Danio rerio
    "NCBITaxon:9606",  # Homo sapiens
    "NCBITaxon:10090",  # Mus musculus
}


# ============================================================================
# Schema Compliance Utilities
# ============================================================================


def _convert_to_ensembl_format(chromosome: str) -> str:
    """Convert chromosome name to ENSEMBL format.

    Per schema 7.1.0:
    - Remove "chr" prefix for GENCODE/UCSC identifiers
    - Convert mitochondrial "M" to "MT"

    :param str chromosome: Chromosome name (e.g., "chr1", "chrM", "1", "MT")
    :return: ENSEMBL-formatted chromosome name
    :rtype: str
    """
    chrom = chromosome

    # Remove "chr" prefix if present
    if chrom.startswith("chr"):
        chrom = chrom[3:]

    # Convert M to MT for mitochondrial chromosome
    if chrom == "M":
        chrom = "MT"

    return chrom


def _remove_ensembl_version(gene_id: str) -> str:
    """Remove version number from Ensembl gene ID.

    Per schema 7.1.0:
    Version numbers MUST be removed from gene_id if it is prefixed with "ENS".
    Example: "ENSG00000186092.7" → "ENSG00000186092"

    :param str gene_id: Gene ID possibly with version number
    :return: Gene ID without version number
    :rtype: str
    """
    # Only remove version for Ensembl IDs (starting with "ENS")
    if gene_id.startswith("ENS") and "." in gene_id:
        return gene_id.split(".")[0]
    return gene_id


# ============================================================================
# Data Processing & Formatting
# ============================================================================


def get_best_matches(csv_file: str) -> pd.DataFrame:
    """Process guidescan2 enumerate output to get the best match per guide.

    - Keep best match per guide_id based on distance
    - Prefer non-NA matches over NA matches
    - Among non-NA matches, prefer lowest distance

    :param str csv_file: Path to guidescan2 enumerate CSV output
    :return: DataFrame with best match per guide_id
    :rtype: pd.DataFrame
    """
    # Load the entire CSV into a DataFrame
    # Keep "NA" as a string, not as pandas NA
    df = pd.read_csv(csv_file, keep_default_na=False, na_values=[""])

    # Create helper column: True for NA matches, False for actual matches
    df["is_no_match"] = df["chromosome"] == "NA"

    # Convert distance to numeric, treating NA as infinity
    df["distance"] = pd.to_numeric(df["distance"], errors="coerce").fillna(float("inf"))

    # Sort by: guide_id, is_no_match (False first), then distance (ascending)
    # This ensures non-NA matches come before NA matches, and lower distances come first
    df_sorted = df.sort_values(by=["id", "is_no_match", "distance"])

    # Group by guide_id and take the first row (which is the best match)
    best_matches = df_sorted.groupby("id", as_index=False).first()

    # Drop the helper column
    best_matches = best_matches.drop(columns=["is_no_match"])

    return best_matches


def split_exact_and_no_matches(best_matches_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Split best matches into exact matches and no matches.

    :param pd.DataFrame best_matches_df: DataFrame with best matches from get_best_matches()
    :return: Tuple of (exact_matches, no_matches) where:

        - exact_matches: distance == 0 AND chromosome != "NA"
        - no_matches: chromosome == "NA"

    :rtype: Tuple[pd.DataFrame, pd.DataFrame]
    """
    # Exact matches: distance == 0 AND chromosome != "NA"
    # Note: distance is already numeric from get_best_matches()
    exact_matches = best_matches_df[(best_matches_df["distance"] == 0) & (best_matches_df["chromosome"] != "NA")].copy()

    # No matches: chromosome == "NA"
    no_matches = best_matches_df[best_matches_df["chromosome"] == "NA"].copy()

    return exact_matches, no_matches


def format_exact_matches(exact_df: pd.DataFrame) -> pd.DataFrame:
    """Format exact matches for bioframe processing.

    Extracts protospacer and PAM from full sequence, converts coordinates
    to BED format (0-based start, 1-based end).

    :param pd.DataFrame exact_df: DataFrame with exact matches
    :return: Formatted DataFrame with columns: id, sequence, pam, chromosome, start, end, sense
    :rtype: pd.DataFrame
    """
    if exact_df.empty:
        return pd.DataFrame(columns=["id", "sequence", "pam", "chromosome", "start", "end", "sense"])

    formatted = []

    for _, row in exact_df.iterrows():
        guide_id = row["id"]
        full_sequence = row["sequence"]
        chromosome = row["chromosome"]
        position = int(row["position"])
        strand = row["strand"]

        # Parse sequence and PAM
        pam_length = len(str(row["pam"]))
        sequence_length = len(full_sequence)
        protospacer_length = sequence_length - pam_length

        protospacer = full_sequence[:protospacer_length]
        pam = full_sequence[protospacer_length:]

        # Convert coordinates to BED format
        # guidescan2 uses 1-based inclusive coordinates
        # BED format uses 0-based start, 1-based end
        if strand == "+":
            # For + strand, position is the start
            start = position - 1  # Convert to 0-based
            end = position + protospacer_length - 1  # BED end coordinate
        else:
            # For - strand, position is the end
            end = position
            start = position - protospacer_length  # 0-based

        formatted.append(
            {
                "id": guide_id,
                "sequence": protospacer,
                "pam": pam,
                "chromosome": chromosome,
                "start": start,
                "end": end,
                "sense": strand,
            }
        )

    return pd.DataFrame(formatted)


def load_gene_coordinates(species: str) -> pd.DataFrame:
    """Load preprocessed gene coordinates for a species.

    :param str species: Species description (e.g., "homo_sapiens", "mus_musculus") or
        NCBITaxon ID (e.g., "NCBITaxon:9606")
    :return: DataFrame with gene coordinates in bioframe format
    :rtype: pd.DataFrame
    :raises FileNotFoundError: If gene coordinates file doesn't exist
    :raises ValueError: If species format is invalid
    """
    # Handle both NCBITaxon format and direct species name
    if species.startswith("NCBITaxon:"):
        # Convert NCBITaxon ID to species name
        try:
            organism = SupportedOrganisms(species)
            # Simply use the enum name in lowercase
            species_name = organism.name.lower()
        except ValueError:
            raise ValueError(f"Invalid NCBITaxon ID: {species}") from None
    else:
        species_name = species

    # Load gene coordinates file
    coord_file = os.path.join(env.GENCODE_DIR, f"gene_coordinates_{species_name}.csv.gz")

    if not os.path.exists(coord_file):
        raise FileNotFoundError(
            f"Gene coordinates file not found: {coord_file}\n"
            f"Please run preprocess_gtf_for_annotation.py first to generate gene coordinates."
        )

    # Load with appropriate column names for bioframe (already uses 'chrom')
    # Specify dtypes for efficiency and type safety
    genes_df = pd.read_csv(coord_file, dtype={"start": "int64", "end": "int64"})

    return genes_df


def normalize_chromosomes(guides_df: pd.DataFrame, genes_df: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Normalize chromosome naming between guides and genes DataFrames.

    Auto-detects chromosome format (chr1 vs 1) and normalizes both
    DataFrames to match.

    :param pd.DataFrame guides_df: DataFrame with guide coordinates
    :param pd.DataFrame genes_df: DataFrame with gene coordinates
    :return: Tuple of (normalized_guides_df, normalized_genes_df)
    :rtype: Tuple[pd.DataFrame, pd.DataFrame]
    """
    if guides_df.empty:
        return guides_df, genes_df

    # Check if guides use "chr" prefix
    first_guide_chr = str(guides_df.iloc[0]["chromosome"])
    guides_use_chr = first_guide_chr.startswith("chr")

    # Check if genes use "chr" prefix
    first_gene_chr = str(genes_df.iloc[0]["chrom"])
    genes_use_chr = first_gene_chr.startswith("chr")

    # Normalize if they differ
    if guides_use_chr and not genes_use_chr:
        # Add "chr" prefix to genes
        genes_df = genes_df.copy()
        genes_df["chrom"] = "chr" + genes_df["chrom"].astype(str)
    elif not guides_use_chr and genes_use_chr:
        # Remove "chr" prefix from genes
        genes_df = genes_df.copy()
        genes_df["chrom"] = genes_df["chrom"].astype(str).str.replace("^chr", "", regex=True)

    return guides_df, genes_df


# ============================================================================
# Bioframe Operations
# ============================================================================


def find_gene_overlaps(guides_df: pd.DataFrame, genes_df: pd.DataFrame) -> pd.DataFrame:
    """Find gene overlaps for guides using bioframe.

    :param pd.DataFrame guides_df: DataFrame with guide coordinates
        (must have: chromosome, start, end)
    :param pd.DataFrame genes_df: DataFrame with gene coordinates
        (must have: chrom, start, end, gene_id, gene_name)
    :return: DataFrame with one row per guide-gene overlap
    :rtype: pd.DataFrame
    """
    if guides_df.empty:
        return pd.DataFrame()

    # Normalize chromosome naming
    guides_df, genes_df = normalize_chromosomes(guides_df, genes_df)

    # Prepare guides for bioframe (rename chromosome to chrom)
    guides_for_overlap = guides_df.copy()
    guides_for_overlap = guides_for_overlap.rename(columns={"chromosome": "chrom"})

    # Ensure proper data types
    guides_for_overlap["start"] = guides_for_overlap["start"].astype("int64")
    guides_for_overlap["end"] = guides_for_overlap["end"].astype("int64")

    # Find overlaps using bioframe
    # how='inner' returns only intervals that overlap (≥1bp)
    overlaps = bf.overlap(
        genes_df, guides_for_overlap, how="inner", cols1=["chrom", "start", "end"], cols2=["chrom", "start", "end"]
    )

    return overlaps


def create_annotated_output(
    overlaps_df: pd.DataFrame, guides_df: pd.DataFrame, no_match_df: pd.DataFrame
) -> pd.DataFrame:
    """Create final annotated output combining overlaps, non-overlapping guides, and no-match guides.

    :param pd.DataFrame overlaps_df: DataFrame with guide-gene overlaps from bioframe
    :param pd.DataFrame guides_df: DataFrame with all formatted guides
    :param pd.DataFrame no_match_df: DataFrame with guides that had no genomic match
    :return: Final annotated output with columns:
        id, sequence, pam, chromosome, start, end, sense, gene_id, gene_name
    :rtype: pd.DataFrame
    """
    result_dfs = []

    # Process overlaps using vectorized operations
    if not overlaps_df.empty:
        # bioframe adds '_' suffix to columns from second dataframe (guides)
        # Keep only the columns we need and rename guide columns
        overlaps_formatted = overlaps_df[
            ["id_", "sequence_", "pam_", "chrom_", "start_", "end_", "sense_", "gene_id", "gene_name"]
        ].rename(
            columns={
                "id_": "id",
                "sequence_": "sequence",
                "pam_": "pam",
                "chrom_": "chromosome",
                "start_": "start",
                "end_": "end",
                "sense_": "sense",
            }
        )
        result_dfs.append(overlaps_formatted)

    # Find guides with no overlaps (guides that were exact matches but didn't overlap any genes)
    if not guides_df.empty:
        overlapping_guide_ids = set(overlaps_df["id_"].unique()) if not overlaps_df.empty else set()
        non_overlapping_guides = guides_df[~guides_df["id"].isin(overlapping_guide_ids)].copy()
        non_overlapping_guides["gene_id"] = "NA"
        non_overlapping_guides["gene_name"] = "NA"
        result_dfs.append(non_overlapping_guides)

    # Add no-match guides
    if not no_match_df.empty:
        no_match_formatted = pd.DataFrame(
            {
                "id": no_match_df["id"],
                "sequence": no_match_df["sequence"],
                "pam": "",  # No PAM info for no-match guides
                "chromosome": "NA",
                "start": "NA",
                "end": "NA",
                "sense": "NA",
                "gene_id": "NA",
                "gene_name": "NA",
            }
        )
        result_dfs.append(no_match_formatted)

    return pd.concat(result_dfs, ignore_index=True) if result_dfs else pd.DataFrame()


# ============================================================================
# Manage GuideScan2 Index Files
# ============================================================================


def _get_guidescan_index_url(organism_id: str) -> str:
    """Get download URL for guidescan2 index (STUB - to be implemented).

    :param str organism_id: NCBITaxon ID (e.g., "NCBITaxon:9606")
    :return: URL to download index tar.gz file
    :rtype: str
    :raises NotImplementedError: This is a stub for future implementation
    """
    # TODO: Map organism IDs to index URLs
    # Will need to host index files and provide download URLs
    raise NotImplementedError(
        "GuideScan2 index download not yet implemented. " "Please provide index path using --index-path option."
    )


def _download_and_extract_index(organism_id: str, cache_dir: str) -> str:
    """Download and extract guidescan2 index tar.gz file.

    :param str organism_id: NCBITaxon ID (e.g., "NCBITaxon:9606")
    :param str cache_dir: Directory to cache downloaded index
    :return: Path to index_prefix (without file extension)
    :rtype: str
    :raises NotImplementedError: If index URL is not available (via _get_guidescan_index_url)
    """
    # Get organism name for file naming
    organism = SupportedOrganisms(organism_id)
    species_name = organism.name.lower()

    # Get download URL
    url = _get_guidescan_index_url(organism_id)

    # Download tar.gz file
    tar_path = os.path.join(cache_dir, f"guidescan_index_{species_name}.tar.gz")
    logger.info(f"Downloading index from {url} to {tar_path}...")

    # Implement download logic
    urllib.request.urlretrieve(url, tar_path)

    # Extract tar.gz safely
    logger.info(f"Extracting {tar_path}...")
    with tarfile.open(tar_path, "r:gz") as tar:
        # Validate all members are safe (no path traversal)
        for member in tar.getmembers():
            if member.name.startswith("/") or ".." in member.name:
                raise ValueError(f"Unsafe tar member: {member.name}")
        tar.extractall(cache_dir)

    # Return index prefix path
    index_prefix = os.path.join(cache_dir, f"guidescan_index_{species_name}")
    logger.info(f"Index extracted to {index_prefix}")

    return index_prefix


def _get_or_download_index(organism_id: str, index_path: Optional[str] = None) -> str:
    """Get cached index or download if not present.

    :param str organism_id: NCBITaxon ID (e.g., "NCBITaxon:9606")
    :param Optional[str] index_path: Explicit path to index. If None, uses cache.
    :return: Path to index_prefix (without file extension)
    :rtype: str
    """
    # If explicit path provided, use it
    if index_path:
        logger.info(f"Using provided index path: {index_path}")
        return index_path

    # Otherwise, check cache or download
    organism = SupportedOrganisms(organism_id)
    species_name = organism.name.lower()

    # Check if index exists in cache
    cache_dir = env.GENCODE_DIR
    index_prefix = os.path.join(cache_dir, f"guidescan_index_{species_name}")

    # Check for any index files with this prefix
    # Guidescan2 indices consist of 3 files with the extensions *.index.{forward,reverse,gs}
    required_extensions = [".index.forward", ".index.reverse", ".index.gs"]
    if all(os.path.exists(f"{index_prefix}{ext}") for ext in required_extensions):
        logger.info(f"Using cached index: {index_prefix}")
        return index_prefix

    # Index not found, attempt download
    logger.info("Index not found in cache, attempting download...")
    return _download_and_extract_index(organism_id, cache_dir)


# ============================================================================
# GuideScan2 Integration
# ============================================================================


def _check_guidescan2_installed() -> None:
    """Check if guidescan2 is installed and available.

    :raises RuntimeError: If guidescan command is not found
    """
    logger.info("Checking for guidescan2 installation...")
    if shutil.which("guidescan") is None:
        raise RuntimeError(
            "guidescan2 not found. Please install guidescan-cli from " "https://github.com/pritykinlab/guidescan-cli"
        )
    logger.info("✓ guidescan2 found")


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


# ============================================================================
# AnnData Extraction & Validation
# ============================================================================


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
            f"Supported organisms: {', '.join(sorted(SUPPORTED_PERTURBATION_ORGANISMS))}"
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
    them to a CSV file in guidescan2 enumerate input format (id, sequence, pam, chromosome, start, end, sense).
    Only id, sequence, and pam are required to run guidescan enumerate. The other columns are empty.

    :param ad.AnnData adata: AnnData object with genetic_perturbations in uns
    :param str output_csv: Path to output CSV file
    """
    genetic_perturbations = adata.uns["genetic_perturbations"]

    # Extract guide sequences
    rows = []
    for guide_id, guide_data in genetic_perturbations.items():
        # Get protospacer and PAM
        protospacer = guide_data.get("protospacer_sequence", "")
        pam_motif = guide_data.get("protospacer_adjacent_motif", "")

        # Extract PAM sequence (remove "3' " prefix if present)
        pam = pam_motif[3:] if pam_motif.startswith("3' ") else pam_motif

        # Combine into full sequence for guidescan2
        full_sequence = protospacer + pam

        rows.append(
            {
                "id": guide_id,
                "sequence": full_sequence,
                "pam": pam,
                "chromosome": "",
                "start": "",
                "end": "",
                "sense": "",
            }
        )

    # Write to CSV
    df = pd.DataFrame(rows)
    df.to_csv(output_csv, index=False)
    logger.info(f"Extracted {len(rows)} guide sequences to {output_csv}")


# ============================================================================
# AnnData Update Helpers
# ============================================================================


def _build_genomic_regions(valid_rows: pd.DataFrame) -> List[str]:
    """Build list of genomic region strings from annotation rows.

    Converts BED format coordinates to schema 7.1.0 format (1-based, fully-closed).

    :param pd.DataFrame valid_rows: DataFrame with columns: chromosome, start, end, sense
    :return: List of genomic region strings in format "chrom:start-end(sense)"
    :rtype: List[str]
    """
    regions = set()
    for _, row in valid_rows.iterrows():
        if row["chromosome"] == "NA":
            continue

        # Convert chromosome to ENSEMBL format
        chrom = _convert_to_ensembl_format(str(row["chromosome"]))

        # Convert from BED (0-based start, 1-based end) to 1-based fully-closed
        start_1based = int(row["start"]) + 1  # Convert 0-based to 1-based
        end_1based = int(row["end"])  # BED end equals 1-based inclusive numerically

        region = f"{chrom}:{start_1based}-{end_1based}({row['sense']})"
        regions.add(region)

    return list(regions)


def _build_target_features(valid_rows: pd.DataFrame) -> Dict[str, str]:
    """Build target_features dict from annotation rows.

    Removes version numbers from Ensembl IDs per schema 7.1.0.

    :param pd.DataFrame valid_rows: DataFrame with columns: gene_id, gene_name
    :return: Dictionary mapping gene_id to gene_name
    :rtype: Dict[str, str]
    """
    features: Dict[str, str] = {}
    for _, row in valid_rows.iterrows():
        if row["gene_id"] != "NA":
            # Remove version number from Ensembl IDs
            gene_id_no_version = _remove_ensembl_version(row["gene_id"])
            features[gene_id_no_version] = row["gene_name"]

    return features


def update_h5ad_with_guide_annotations(adata: ad.AnnData, annotations_df: pd.DataFrame) -> ad.AnnData:
    """Update h5ad genetic_perturbations with genomic annotations from guidescan2/bioframe.

    Takes the output from annotate_guides_from_guidescan_csv and adds genomic location
    and gene overlap information to adata.uns['genetic_perturbations']. Only updates
    target_genomic_regions and target_features fields (does NOT modify protospacer_sequence
    or protospacer_adjacent_motif which should already be present in the h5ad).

    :param ad.AnnData adata: AnnData object with genetic_perturbations in uns
    :param pd.DataFrame annotations_df: DataFrame from annotate_guides_from_guidescan_csv with
        required columns: id, chromosome, start, end, sense, gene_id, gene_name
    :return: Modified AnnData object
    :rtype: ad.AnnData
    :raises KeyError: If genetic_perturbations is not in adata.uns
    """
    # Check that genetic_perturbations exists
    if "genetic_perturbations" not in adata.uns:
        raise KeyError("adata.uns does not contain 'genetic_perturbations'")

    genetic_perturbations = adata.uns["genetic_perturbations"]

    # Filter annotations to only guides that exist in h5ad
    h5ad_guide_ids = set(genetic_perturbations.keys())
    extra_guides = set(annotations_df["id"].unique()) - h5ad_guide_ids

    if extra_guides:
        logger.warning(
            f"Skipping {len(extra_guides)} guide(s) from annotations that are not in h5ad: " f"{sorted(extra_guides)}"
        )

    # Filter to valid guides with gene overlaps (vectorized)
    valid_annotations = annotations_df[
        (annotations_df["id"].isin(h5ad_guide_ids)) & (annotations_df["gene_id"] != "NA")
    ].copy()

    # Process each guide using groupby
    for guide_id, guide_group in valid_annotations.groupby("id"):
        # Build genomic regions and target features
        genomic_regions = _build_genomic_regions(guide_group)
        target_features = _build_target_features(guide_group)

        # Update the h5ad entry
        genetic_perturbations[guide_id]["target_genomic_regions"] = genomic_regions
        genetic_perturbations[guide_id]["target_features"] = target_features

        logger.info(f"Updated guide {guide_id}: {len(genomic_regions)} region(s), {len(target_features)} gene(s)")

    # Update the uns with modified genetic_perturbations
    adata.uns["genetic_perturbations"] = genetic_perturbations

    return adata


# ============================================================================
# High-Level Pipeline Functions (Public API)
# ============================================================================


def annotate_guides_from_guidescan_csv(
    guidescan_csv: str, species: str, output_csv: Optional[str] = None
) -> pd.DataFrame:
    """Main function to annotate guides from guidescan2 enumerate output.

    :param str guidescan_csv: Path to guidescan2 enumerate CSV output
    :param str species: Species NCBITaxon ID (e.g., "NCBITaxon:9606") or
        species name (e.g., "homo_sapiens")
    :param Optional[str] output_csv: Path to save annotated output CSV. If None, doesn't save.
    :return: Annotated guides DataFrame
    :rtype: pd.DataFrame
    """
    # Step 1: Get best matches
    logger.info("Processing best matches...")
    best_matches = get_best_matches(guidescan_csv)
    logger.info(f"Found {len(best_matches)} unique guides")

    # Step 2: Split exact vs no matches
    logger.info("Splitting exact and no matches...")
    exact_matches, no_matches = split_exact_and_no_matches(best_matches)
    logger.info(f"  Exact matches: {len(exact_matches)}")
    logger.info(f"  No matches: {len(no_matches)}")

    # Step 3: Format exact matches
    logger.info("Formatting exact matches...")
    formatted_guides = format_exact_matches(exact_matches)

    # Step 4: Load gene coordinates
    logger.info(f"Loading gene coordinates for {species}...")
    genes_df = load_gene_coordinates(species)
    logger.info(f"Loaded {len(genes_df)} genes")

    # Step 5: Find overlaps
    logger.info("Finding gene overlaps...")
    overlaps = find_gene_overlaps(formatted_guides, genes_df)
    logger.info(f"Found {len(overlaps)} guide-gene overlaps")

    # Step 6: Create final output
    logger.info("Creating annotated output...")
    annotated = create_annotated_output(overlaps, formatted_guides, no_matches)
    logger.info(f"Final output: {len(annotated)} rows")

    # Save if output path provided
    if output_csv:
        logger.info(f"Saving to {output_csv}...")
        annotated.to_csv(output_csv, index=False)
        logger.info("Done!")

    return annotated


def annotate_perturbations_in_h5ad(adata: ad.AnnData, index_path: Optional[str] = None) -> ad.AnnData:
    """Annotate genetic perturbations in an AnnData object.

    This function orchestrates the annotation pipeline:
    1. Check if genetic_perturbations exists (skip if not)
    2. Check guidescan2 is installed
    3. Detect organism from adata.uns
    4. Extract perturbations to temp CSV
    5. Get/download guidescan2 index for organism
    6. Run guidescan enumerate (--mismatches 0)
    7. Run bioframe annotation
    8. Update adata with annotations
    9. Return modified adata

    :param ad.AnnData adata: AnnData object (optionally with genetic_perturbations in uns)
    :param Optional[str] index_path: Path to guidescan2 index. If None, uses default cache.
    :return: Modified AnnData object with annotations (or unmodified if no perturbations)
    :rtype: ad.AnnData
    :raises RuntimeError: If guidescan2 is not installed
    """
    # Check if genetic_perturbations exists
    if "genetic_perturbations" not in adata.uns:
        logger.info("No genetic_perturbations found in adata.uns, skipping annotation")
        return adata

    # Check guidescan2 is installed
    _check_guidescan2_installed()

    # Detect organism
    logger.info("Detecting organism...")
    organism_id = _get_organism_from_h5ad(adata)

    # Use temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmpdir:
        # Extract perturbations to temp CSV
        logger.info("Extracting perturbations to CSV...")
        guides_csv = os.path.join(tmpdir, "guides_input.csv")
        _extract_perturbations_to_guidescan_csv(adata, guides_csv)

        # Get/download guidescan2 index
        logger.info("Getting guidescan2 index...")
        try:
            index_prefix = _get_or_download_index(organism_id, index_path)
        except NotImplementedError as e:
            logger.error(str(e))
            raise RuntimeError(
                f"GuideScan2 index for {organism_id} not found. "
                "Please provide --index-path or implement index download."
            ) from e

        # Run guidescan enumerate
        logger.info("Running guidescan enumerate...")
        enumerate_csv = os.path.join(tmpdir, "guidescan_enumerate_output.csv")
        _run_guidescan_enumerate(guides_csv, index_prefix, enumerate_csv)

        # Run bioframe annotation
        logger.info("Annotating with bioframe...")
        annotations_df = annotate_guides_from_guidescan_csv(enumerate_csv, organism_id)

        # Update adata with annotations
        logger.info("Updating AnnData with annotations...")
        adata = update_h5ad_with_guide_annotations(adata, annotations_df)

    logger.info("✓ Annotation complete!")
    return adata
