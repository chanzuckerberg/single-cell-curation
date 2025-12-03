"""Annotate guidescan2 gRNA output with gene overlaps using bioframe.

This module processes guidescan2 enumerate output and annotates gRNAs with
information about which genes they overlap using genomic interval analysis.
"""

import logging
import os
from typing import Optional, Tuple

import anndata as ad
import bioframe as bf
import pandas as pd

from . import env
from .gencode import SupportedOrganisms

logger = logging.getLogger(__name__)


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
        # PAM is typically last 3 nucleotides (e.g., NGG for Cas9)
        pam_length = 3
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
            end = position + protospacer_length - 1  # Keep as 1-based
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
    annotated_rows = []

    if not overlaps_df.empty:
        # Process overlaps
        # bioframe adds '_' suffix to columns from second dataframe
        for _, row in overlaps_df.iterrows():
            annotated_rows.append(
                {
                    "id": row["id_"],
                    "sequence": row["sequence_"],
                    "pam": row["pam_"],
                    "chromosome": row["chrom_"],
                    "start": row["start_"],
                    "end": row["end_"],
                    "sense": row["sense_"],
                    "gene_id": row["gene_id"],
                    "gene_name": row["gene_name"],
                }
            )

    # Find guides with no overlaps (guides that were exact matches but didn't overlap any genes)
    if not guides_df.empty:
        overlapping_guide_ids = set(overlaps_df["id_"].unique()) if not overlaps_df.empty else set()

        non_overlapping_guides = guides_df[~guides_df["id"].isin(overlapping_guide_ids)]

        for _, row in non_overlapping_guides.iterrows():
            annotated_rows.append(
                {
                    "id": row["id"],
                    "sequence": row["sequence"],
                    "pam": row["pam"],
                    "chromosome": row["chromosome"],
                    "start": row["start"],
                    "end": row["end"],
                    "sense": row["sense"],
                    "gene_id": "NA",
                    "gene_name": "NA",
                }
            )

    # Add no-match guides
    if not no_match_df.empty:
        for _, row in no_match_df.iterrows():
            annotated_rows.append(
                {
                    "id": row["id"],
                    "sequence": row["sequence"],
                    "pam": "",  # No PAM info for no-match guides
                    "chromosome": "NA",
                    "start": "NA",
                    "end": "NA",
                    "sense": "NA",
                    "gene_id": "NA",
                    "gene_name": "NA",
                }
            )

    return pd.DataFrame(annotated_rows)


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
    :raises ValueError: If any guide IDs from annotations_df are missing from genetic_perturbations
    :raises KeyError: If genetic_perturbations is not in adata.uns
    """
    # Check that genetic_perturbations exists
    if "genetic_perturbations" not in adata.uns:
        raise KeyError("adata.uns does not contain 'genetic_perturbations'")

    genetic_perturbations = adata.uns["genetic_perturbations"]

    # Validate that all guides in annotations exist in h5ad
    annotation_guide_ids = set(annotations_df["id"].unique())
    h5ad_guide_ids = set(genetic_perturbations.keys())
    missing_guides = annotation_guide_ids - h5ad_guide_ids

    if missing_guides:
        raise ValueError(
            f"The following guide IDs from annotations are missing from "
            f"adata.uns['genetic_perturbations']: {sorted(missing_guides)}"
        )

    # Process each guide
    for guide_id in annotation_guide_ids:
        # Get all rows for this guide
        guide_rows = annotations_df[annotations_df["id"] == guide_id]

        # Skip guides with no gene overlaps (all gene_id are "NA")
        if (guide_rows["gene_id"] == "NA").all():
            logger.debug(f"Skipping guide {guide_id}: no gene overlaps found")
            continue

        # Filter out NA rows for processing
        valid_rows = guide_rows[guide_rows["gene_id"] != "NA"]

        # Build target_genomic_regions: unique list of "chrom:start-end(sense)"
        # Schema requires 1-based, fully-closed coordinates (not BED format)
        genomic_regions = set()
        for _, row in valid_rows.iterrows():
            # Skip rows with NA coordinates
            if row["chromosome"] == "NA":
                continue

            # Convert chromosome to ENSEMBL format
            chrom = _convert_to_ensembl_format(str(row["chromosome"]))

            # Convert from BED (0-based start, 1-based end) to 1-based fully-closed
            start_1based = int(row["start"]) + 1  # Convert 0-based to 1-based
            end_1based = int(row["end"])  # Already 1-based in BED format

            region = f"{chrom}:{start_1based}-{end_1based}({row['sense']})"
            genomic_regions.add(region)

        # Build target_features: dict of {gene_id: gene_name}
        # Schema requires removing version numbers from Ensembl IDs
        target_features = {}
        for _, row in valid_rows.iterrows():
            gene_id = row["gene_id"]
            gene_name = row["gene_name"]
            if gene_id != "NA":
                # Remove version number from Ensembl IDs (e.g., "ENSG00000186092.7" → "ENSG00000186092")
                gene_id_no_version = _remove_ensembl_version(gene_id)
                target_features[gene_id_no_version] = gene_name

        # Update the h5ad entry
        genetic_perturbations[guide_id]["target_genomic_regions"] = list(genomic_regions)
        genetic_perturbations[guide_id]["target_features"] = target_features

        logger.info(f"Updated guide {guide_id}: {len(genomic_regions)} region(s), " f"{len(target_features)} gene(s)")

    # Update the uns with modified genetic_perturbations
    adata.uns["genetic_perturbations"] = genetic_perturbations

    return adata
