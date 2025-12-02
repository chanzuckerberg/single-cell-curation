#!/usr/bin/env python3
"""
Extract genetic perturbation data from uns['genetic_perturbations'] and write to CSV.

This script reads an h5ad file and extracts genetic perturbation information
from uns['genetic_perturbations'] into a CSV file with columns:
id, sequence, pam, chromosome, start, end, sense
"""

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import anndata as ad


def parse_genomic_region(region_str: str) -> Optional[Tuple[str, str, str, str]]:
    """
    Parse a genomic region string in the format "CHROMOSOME:START-END(STRAND)".

    Parameters
    ----------
    region_str : str
        Genomic region string, e.g., "1:1000-1018(+)"

    Returns
    -------
    tuple or None
        (chromosome, start, end, sense) or None if parsing fails
    """
    # Pattern: CHROMOSOME:START-END(STRAND)
    # Examples: "1:1000-1018(+)", "MT:500-520(-)"
    pattern = r"^([^:]+):(\d+)-(\d+)\(([+-])\)$"
    match = re.match(pattern, region_str)
    if match:
        chromosome = match.group(1)
        start = match.group(2)
        end = match.group(3)
        sense = match.group(4)
        return (chromosome, start, end, sense)
    return None


def extract_pam_from_motif(motif: str) -> str:
    """
    Extract PAM sequence from protospacer_adjacent_motif.

    The format is "3' MOTIF", e.g., "3' NGG" -> "NGG"

    Parameters
    ----------
    motif : str
        protospacer_adjacent_motif value

    Returns
    -------
    str
        PAM sequence without the "3' " prefix
    """
    # Remove "3' " prefix if present
    if motif.startswith("3' "):
        return motif[3:]
    return motif


def extract_perturbations_to_csv(
    input_file: str,
    output_file: str,
) -> None:
    """
    Extract genetic perturbation data from h5ad file and write to CSV.

    Parameters
    ----------
    input_file : str
        Path to input h5ad file
    output_file : str
        Path to output CSV file
    """
    print(f"Reading {input_file}...")
    try:
        adata = ad.read_h5ad(input_file)
    except Exception as e:
        print(f"Error reading file: {e}", file=sys.stderr)
        sys.exit(1)

    # Check if genetic_perturbations exists
    if "genetic_perturbations" not in adata.uns:
        print("Error: 'genetic_perturbations' not found in uns", file=sys.stderr)
        sys.exit(1)

    genetic_perturbations: Dict = adata.uns["genetic_perturbations"]
    print(f"Found {len(genetic_perturbations)} genetic perturbation(s)")

    # Prepare rows for CSV
    rows: List[List[str]] = []

    for pert_id, pert_data in genetic_perturbations.items():
        # Extract basic fields
        sequence = pert_data.get("protospacer_sequence", "")
        motif = pert_data.get("protospacer_adjacent_motif", "")
        pam = extract_pam_from_motif(motif)

        # Extract genomic regions
        target_regions = pert_data.get("target_genomic_regions", [])
        if not target_regions:
            # If no target regions, create a row with empty chromosome/position/sense
            rows.append([pert_id, sequence, pam, "", "", "", ""])
        else:
            # Create a row for each target genomic region
            for region_str in target_regions:
                parsed = parse_genomic_region(region_str)
                if parsed:
                    chromosome, start, end, sense = parsed
                    rows.append([pert_id, sequence, pam, chromosome, start, end, sense])
                else:
                    # If parsing fails, include the raw region string
                    print(f"Warning: Could not parse genomic region '{region_str}' for {pert_id}", file=sys.stderr)
                    rows.append([pert_id, sequence, pam, "", "", "", region_str])

    # Write to CSV
    print(f"Writing {len(rows)} row(s) to {output_file}...")
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        # Write header
        writer.writerow(["id", "sequence", "pam", "chromosome", "start", "end", "sense"])
        # Write data rows
        writer.writerows(rows)

    print(f"Done! Wrote {len(rows)} row(s) to {output_file}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract genetic perturbation data from h5ad file to CSV")
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to input h5ad file",
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to output CSV file",
    )

    args = parser.parse_args()

    # Check if input file exists
    if not Path(args.input_file).exists():
        print(f"Error: Input file '{args.input_file}' does not exist", file=sys.stderr)
        sys.exit(1)

    extract_perturbations_to_csv(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
