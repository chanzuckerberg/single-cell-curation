#!/usr/bin/env python3
"""
Create a smaller test h5ad dataset from curated.h5ad that preserves
the genetic perturbation fields required for schema 7.1.0 testing.

This script:
1. Reads curated.h5ad in backed mode to avoid loading everything into memory
2. Creates a subset with fewer cells and genes
3. Preserves all genetic perturbation fields (obs and uns)
4. Ensures the subset maintains schema 7.1.0 compliance
"""

import argparse
import sys
from pathlib import Path
from typing import List

import anndata as ad
import numpy as np
from scipy import sparse


def create_small_test_dataset(
    input_file: str,
    output_file: str,
    n_cells: int = 500,
    n_genes: int = 500,
    n_genetic_perturbations: int = 10,
    seed: int = 42,
) -> None:
    """
    Create a smaller test dataset from a large h5ad file.

    Parameters
    ----------
    input_file : str
        Path to the input h5ad file
    output_file : str
        Path to the output h5ad file
    n_cells : int
        Number of cells to include in the subset (default: 500)
    n_genes : int
        Number of genes to include in the subset (default: 2000)
    seed : int
        Random seed for reproducibility (default: 42)
    """
    print(f"Reading {input_file}...")

    # Read in backed mode to avoid loading everything into memory
    try:
        adata = ad.read_h5ad(input_file, backed="r")
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        print("Trying to read without backed mode...")
        adata = ad.read_h5ad(input_file)

    print(f"Original dataset: {adata.n_obs} cells, {adata.n_vars} genes")

    # Set random seed for reproducibility
    np.random.seed(seed)

    # Determine which cells to keep
    # If we have genetic perturbation data, try to preserve diversity
    if "genetic_perturbation_id" in adata.obs.columns:
        # Sample cells to preserve different perturbation IDs
        perturbation_ids = adata.obs["genetic_perturbation_id"].unique()
        print(f"Found {len(perturbation_ids)} unique genetic perturbation IDs")

        # Sample cells proportionally from each perturbation ID
        cells_to_keep: List[str] = []
        for pert_id in perturbation_ids:
            pert_cells = adata.obs.index[adata.obs["genetic_perturbation_id"] == pert_id]
            n_pert_cells = min(len(pert_cells), max(1, n_cells // len(perturbation_ids)))
            sampled = np.random.choice(pert_cells, size=n_pert_cells, replace=False)
            cells_to_keep.extend(sampled)

        # If we need more cells, randomly sample from the rest
        if len(cells_to_keep) < n_cells:
            remaining_cells = adata.obs.index[~adata.obs.index.isin(cells_to_keep)]
            n_additional = n_cells - len(cells_to_keep)
            if len(remaining_cells) > 0:
                additional = np.random.choice(
                    remaining_cells, size=min(n_additional, len(remaining_cells)), replace=False
                )
                cells_to_keep.extend(additional)

        cells_to_keep = cells_to_keep[:n_cells]
    else:
        # No genetic perturbation data, just random sample
        cells_to_keep = np.random.choice(adata.obs.index, size=min(n_cells, adata.n_obs), replace=False)

    # Determine which genes to keep
    # Keep genes that have non-zero expression in at least some cells
    if adata.n_vars > n_genes:
        print("Calculating gene expression to select top genes...")
        # For backed mode, we need to handle the lazy dataset differently
        # Try to calculate gene sums - if it fails, we're in backed mode
        try:
            # Check if we can call sum() directly
            if sparse.issparse(adata.X):  # noqa
                # Try to get sum - this will fail for backed datasets
                gene_counts = np.array(adata.X.sum(axis=0)).flatten()
            else:
                # Dense array
                gene_counts = adata.X.sum(axis=0)

            # Sort by expression and take top genes
            top_genes_idx = np.argsort(gene_counts)[-n_genes:]
            genes_to_keep = adata.var.index[top_genes_idx]
            print(f"  Selected top {n_genes} genes by expression")
        except (AttributeError, TypeError, ValueError):
            # Fallback: random sample if we can't calculate sums (backed mode)
            print("  Using random sampling strategy for gene selection (backed mode)...")
            genes_to_keep = np.random.choice(adata.var.index, size=min(n_genes, adata.n_vars), replace=False)
    else:
        genes_to_keep = adata.var.index

    print(f"Creating subset: {len(cells_to_keep)} cells, {len(genes_to_keep)} genes")

    # Create the subset
    # First, convert to memory if in backed mode, then subset
    print("Creating subset (this may take a moment for large files)...")
    if hasattr(adata, "filename") and adata.filename is not None:
        # It's in backed mode, convert to memory first
        print("  Converting from backed mode to memory...")
        adata = adata.to_memory()

    # Now subset by cells and genes
    adata_subset = adata[cells_to_keep, genes_to_keep].copy()

    # Ensure we have the genetic perturbation fields in obs
    genetic_pert_fields = ["genetic_perturbation_id", "genetic_perturbation_strategy"]
    for field in genetic_pert_fields:
        if field in adata.obs.columns and field not in adata_subset.obs.columns:
            adata_subset.obs[field] = adata.obs.loc[cells_to_keep, field]

    # Preserve genetic_perturbations in uns
    if "genetic_perturbations" in adata.uns:
        print("Preserving genetic_perturbations in uns...")
        adata_subset.uns["genetic_perturbations"] = adata.uns["genetic_perturbations"].copy()

        # Filter target_features to only include genes that are in our subset
        for pert_id, pert_data in adata_subset.uns["genetic_perturbations"].items():
            if "target_features" in pert_data:
                # Keep only features that are in our gene subset
                target_features = pert_data["target_features"]
                filtered_features = {
                    feat_id: feat_name for feat_id, feat_name in target_features.items() if feat_id in genes_to_keep
                }
                if filtered_features:
                    adata_subset.uns["genetic_perturbations"][pert_id]["target_features"] = filtered_features
                else:
                    # If no features remain, remove target_features key
                    pert_data.pop("target_features", None)
        # Only keep at most n_genetic_perturbations if specified. Select randomly.
        if n_genetic_perturbations > 0:
            pert_ids = list(adata_subset.uns["genetic_perturbations"].keys())
            if len(pert_ids) > n_genetic_perturbations:
                selected_pert_ids = np.random.choice(pert_ids, size=n_genetic_perturbations, replace=False)
                adata_subset.uns["genetic_perturbations"] = {
                    pid: adata_subset.uns["genetic_perturbations"][pid] for pid in selected_pert_ids
                }

    # Ensure raw.var exists if raw exists
    if adata.raw is not None:
        print("Preserving raw data...")
        # Subset raw by the same genes
        raw_subset = adata.raw[:, genes_to_keep].copy()
        adata_subset.raw = raw_subset

    # Ensure obsm embeddings are preserved
    # Note: AnnData automatically subsets obsm when we subset the main object
    if adata.obsm:
        print(f"Preserving {len(adata.obsm)} embedding(s)...")
        # Verify embeddings were subset correctly
        for key in adata_subset.obsm:
            if adata_subset.obsm[key].shape[0] != len(cells_to_keep):
                print(f"  WARNING: {key} shape mismatch ({adata_subset.obsm[key].shape[0]} vs {len(cells_to_keep)})")

    # Ensure schema version is set to 7.1.0
    adata_subset.uns["schema_version"] = "7.1.0"
    adata_subset.uns["schema_reference"] = (
        "https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/7.1.0/schema.md"
    )

    # Verify genetic perturbation fields are present
    print("\nVerifying genetic perturbation fields...")
    has_genetic_pert = False
    if "genetic_perturbation_id" in adata_subset.obs.columns:
        n_with_pert = (adata_subset.obs["genetic_perturbation_id"] != "na").sum()
        n_total = len(adata_subset.obs)
        print(f"  ✓ genetic_perturbation_id: {n_with_pert}/{n_total} cells have perturbation IDs")
        has_genetic_pert = True
    else:
        print("  ⚠ genetic_perturbation_id not found in obs")

    if "genetic_perturbation_strategy" in adata_subset.obs.columns:
        strategies = adata_subset.obs["genetic_perturbation_strategy"].unique()
        print(f"  ✓ genetic_perturbation_strategy: {list(strategies)}")
    else:
        print("  ⚠ genetic_perturbation_strategy not found in obs")

    if "genetic_perturbations" in adata_subset.uns:
        n_perturbations = len(adata_subset.uns["genetic_perturbations"])
        print(f"  ✓ genetic_perturbations: {n_perturbations} perturbation(s) defined in uns")

        # Show details of each perturbation
        for pert_id, pert_data in adata_subset.uns["genetic_perturbations"].items():
            role = pert_data.get("role", "unknown")
            n_targets = len(pert_data.get("target_features", {}))
            print(f"    - {pert_id}: role={role}, {n_targets} target feature(s)")
    else:
        print("  ⚠ genetic_perturbations not found in uns")
        if has_genetic_pert:
            print("    WARNING: obs has genetic_perturbation_id but uns lacks genetic_perturbations!")

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Input:  {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    print(f"Output: {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")
    print(f"Reduction: {adata.n_obs/adata_subset.n_obs:.1f}x cells, {adata.n_vars/adata_subset.n_vars:.1f}x genes")
    print(f"Schema version: {adata_subset.uns.get('schema_version', 'not set')}")
    print(f"Genetic perturbation fields: {'✓ Present' if has_genetic_pert else '✗ Missing'}")
    print("=" * 60)

    # Write the subset
    print(f"\nWriting subset to {output_file}...")
    adata_subset.write_h5ad(output_file, compression="gzip")

    # Print file size
    output_size = Path(output_file).stat().st_size / (1024 * 1024)  # MB
    print(f"Output file size: {output_size:.2f} MB")

    print("\nDone!")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create a smaller test h5ad dataset preserving genetic perturbation fields"
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to input h5ad file (e.g., curated.h5ad)",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to output h5ad file",
    )
    parser.add_argument(
        "--n-cells",
        type=int,
        default=500,
        help="Number of cells to include (default: 500)",
    )
    parser.add_argument(
        "--n-genes",
        type=int,
        default=2000,
        help="Number of genes to include (default: 2000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)",
    )

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    create_small_test_dataset(
        args.input_file,
        args.output_file,
        n_cells=args.n_cells,
        n_genes=args.n_genes,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
