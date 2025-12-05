#!/usr/bin/env python3
"""Preprocess GTF files to create gene coordinate BED files for guidescan2 annotation.

This script downloads GTF files from URLs specified in gene_info.yml and extracts
gene coordinates into BED-format CSV files for fast lookup during gRNA annotation.

By default, only processes species with perturbation_support: true in gene_info.yml
(human, mouse, zebrafish). Use --species to process a specific species regardless
of perturbation support.

Output format: chromosome, start, end, gene_id, gene_name, strand
Coordinates are in BED format (0-based start, 1-based end)

Usage::

    python preprocess_gtf_for_annotation.py [--species SPECIES]
"""

import argparse
import gzip
import os
import sys
import urllib.request
from typing import Dict, List, Tuple

import yaml

# Add cellxgene_schema to path
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), "../cellxgene_schema"))
import env
import gtf_tools


def download_gtf(url: str, output_path: str) -> None:
    """Download a GTF file from a URL.

    :param str url: URL to download from
    :param str output_path: Path to save the downloaded file
    :rtype: None
    """
    if os.path.exists(output_path):
        print(f"GTF file already exists at {output_path}, skipping download")
        return

    print(f"Downloading {url}...")
    urllib.request.urlretrieve(url, output_path)
    print(f"Downloaded to {output_path}")


def parse_gtf_for_gene_coordinates(gtf_path: str) -> List[Tuple[str, int, int, str, str, str]]:
    """Parse a GTF file to extract gene coordinates.

    :param str gtf_path: Path to gzipped GTF file
    :return: List of tuples, each containing (chromosome, start, end, gene_id, gene_name, strand).
        Coordinates are in BED format (0-based start, 1-based end)
    :rtype: List[Tuple[str, int, int, str, str, str]]
    """
    gene_coordinates = []

    with gzip.open(gtf_path, "rb") as gtf:
        for byte_line in gtf:
            line = byte_line.decode("utf-8")

            # Skip comment lines
            if line[0] == "#":
                continue

            line_parts = line.rstrip().split("\t")

            # Only process gene lines
            if line_parts[2] != "gene":
                continue

            # Extract basic fields
            chromosome = line_parts[0]
            start = int(line_parts[3]) - 1  # Convert to 0-based
            end = int(line_parts[4])  # Keep 1-based (BED format)
            strand = line_parts[6]

            # Extract features from column 9
            features = gtf_tools._get_features(line_parts)

            # Get gene_id and gene_name
            gene_id = features.get("gene_id", "")
            gene_name = features.get("gene_name", gene_id)  # Default to gene_id if name not present

            # Skip PAR_Y genes
            if gene_id.endswith("PAR_Y"):
                continue

            gene_coordinates.append((chromosome, start, end, gene_id, gene_name, strand))

    return gene_coordinates


def write_gene_coordinates(coordinates: List[Tuple], output_path: str) -> None:
    """Write gene coordinates to a gzipped CSV file in bioframe-compatible format.

    :param List[Tuple] coordinates: Gene coordinates to write
    :param str output_path: Path to output file
    :rtype: None
    """
    print(f"Writing {len(coordinates)} gene coordinates to {output_path}...")

    # Sort by chromosome and start position
    coordinates.sort(key=lambda x: (x[0], x[1]))

    with gzip.open(output_path, "wt") as f:
        # Write header with bioframe-compatible column name (chrom instead of chromosome)
        f.write("chrom,start,end,gene_id,gene_name,strand\n")

        # Write data
        for chrom, start, end, gene_id, gene_name, strand in coordinates:
            f.write(f"{chrom},{start},{end},{gene_id},{gene_name},{strand}\n")

    print(f"Wrote {len(coordinates)} gene coordinates")


def process_species(species_key: str, species_info: Dict, gencode_dir: str) -> None:
    """Process a single species: download GTF and create gene coordinates file.

    :param str species_key: Species key from gene_info.yml
    :param Dict species_info: Species information including URL and description
    :param str gencode_dir: Directory to store output files
    :rtype: None
    """
    description = species_info.get("description")
    url = species_info.get("url")

    if not description or not url:
        print(f"Skipping {species_key}: missing description or URL")
        return

    print(f"\n{'='*60}")
    print(f"Processing {species_key} ({description})")
    print(f"{'='*60}")

    # Format URL if it contains version placeholder
    version = species_info.get("version")
    if version:
        url = url.format(version=version)

    try:
        # Download GTF file to gencode_dir (will skip if already exists)
        gtf_filename = f"{description}.gtf.gz"
        gtf_path = os.path.join(gencode_dir, gtf_filename)
        download_gtf(url, gtf_path)

        # Parse GTF to extract gene coordinates
        print("Parsing GTF file...")
        coordinates = parse_gtf_for_gene_coordinates(gtf_path)

        # Write gene coordinates to gencode_dir
        output_filename = f"gene_coordinates_{description}.csv.gz"
        output_path = os.path.join(gencode_dir, output_filename)
        write_gene_coordinates(coordinates, output_path)

        print(f"✓ Successfully processed {species_key}")

    except Exception as e:
        print(f"✗ Error processing {species_key}: {e}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preprocess GTF files to create gene coordinate BED files for guidescan2 annotation"
    )
    parser.add_argument(
        "--species",
        type=str,
        default=None,
        help="Process only this species (e.g., 'human', 'mouse'). If not specified, processes all species.",
    )

    args = parser.parse_args()

    # Load gene info
    print(f"Loading gene info from {env.GENE_INFO_YAML}")
    with open(env.GENE_INFO_YAML, "r") as f:
        gene_info = yaml.safe_load(f)

    gencode_dir = env.GENCODE_DIR

    if args.species:
        # Process single species
        if args.species not in gene_info:
            print(f"Error: Species '{args.species}' not found in gene_info.yml")
            print(f"Available species: {', '.join(gene_info.keys())}")
            sys.exit(1)

        process_species(args.species, gene_info[args.species], gencode_dir)
    else:
        # Process only species with perturbation support
        perturbation_species = {key: info for key, info in gene_info.items() if info.get("perturbation_support", False)}
        print(f"Processing {len(perturbation_species)} species with perturbation support...")
        print(f"Species: {', '.join(perturbation_species.keys())}")
        for species_key, species_info in perturbation_species.items():
            process_species(species_key, species_info, gencode_dir)

    print("\n" + "=" * 60)
    print("Preprocessing complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
