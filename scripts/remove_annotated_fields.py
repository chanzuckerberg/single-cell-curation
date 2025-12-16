#!/usr/bin/env python3
"""
Script to remove annotated fields from genetic_perturbations in an h5ad file.

Removes the following fields from uns['genetic_perturbations'][id]:
- 'target_features'
- 'target_genomic_regions'

Usage:
    python scripts/remove_annotated_fields.py input.h5ad output.h5ad
"""

import argparse
import sys
from pathlib import Path

import anndata as ad


def remove_annotated_fields(input_file: str, output_file: str) -> None:
    """
    Remove target_features and target_genomic_regions from genetic_perturbations.

    Args:
        input_file: Path to input h5ad file
        output_file: Path to output h5ad file
    """
    print(f"Loading h5ad file: {input_file}")
    adata = ad.read_h5ad(input_file)

    # Check if genetic_perturbations exists in uns
    if "genetic_perturbations" not in adata.uns:
        print("Warning: 'genetic_perturbations' not found in uns. Nothing to remove.")
        adata.write_h5ad(output_file, compression="gzip")
        return

    genetic_perturbations = adata.uns["genetic_perturbations"]
    removed_count = 0

    # Iterate through all perturbation IDs
    for pert_id in genetic_perturbations:
        pert_data = genetic_perturbations[pert_id]

        # Remove target_features if present
        if "target_features" in pert_data:
            del pert_data["target_features"]
            removed_count += 1
            print(f"Removed 'target_features' from perturbation ID: {pert_id}")

        # Remove target_genomic_regions if present
        if "target_genomic_regions" in pert_data:
            del pert_data["target_genomic_regions"]
            removed_count += 1
            print(f"Removed 'target_genomic_regions' from perturbation ID: {pert_id}")

    if removed_count == 0:
        print("No annotated fields found to remove.")
    else:
        print(f"Removed {removed_count} field(s) total.")

    print(f"Writing modified h5ad file to: {output_file}")
    adata.write_h5ad(output_file, compression="gzip")
    print("Done.")


def main() -> None:
    """Main entry point for the script."""
    parser = argparse.ArgumentParser(
        description="Remove annotated fields (target_features and target_genomic_regions) "
        "from genetic_perturbations in an h5ad file."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to input h5ad file",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to output h5ad file",
    )

    args = parser.parse_args()

    # Validate input file exists
    input_path = Path(args.input_file)
    if not input_path.exists():
        print(f"Error: Input file does not exist: {args.input_file}", file=sys.stderr)
        sys.exit(1)

    # Create output directory if it doesn't exist
    output_path = Path(args.output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        remove_annotated_fields(args.input_file, args.output_file)
    except Exception as e:
        print(f"Error processing file: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
