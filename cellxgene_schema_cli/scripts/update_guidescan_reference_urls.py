#!/usr/bin/env python3
"""Update reference_files.yml with new guidescan2 index S3 URLs.

Derives the expected S3 archive URL for each species in gene_assembly.yml,
mirroring the naming logic in data-pipelines run_guidescan_index.py:

  1. Extract filename from fasta_url  (e.g. GRCh38.primary_assembly.genome.fa.gz)
  2. Strip the .gz extension          (e.g. GRCh38.primary_assembly.genome.fa)
  3. Append .tar.gz                   (e.g. GRCh38.primary_assembly.genome.fa.tar.gz)
  4. Prepend the S3 HTTPS base URL

The script is idempotent: if the URL in reference_files.yml already matches
what is derived from gene_assembly.yml, the file is left unchanged.

Exit codes:
  0 — success (with or without changes)
  1 — unexpected error
"""

import sys
from pathlib import Path

import yaml

GENE_ASSEMBLY_PATH = Path(__file__).parents[1] / "cellxgene_schema" / "gencode_files" / "gene_assembly.yml"
REFERENCE_FILES_PATH = Path(__file__).parents[1] / "cellxgene_schema" / "gencode_files" / "reference_files.yml"

S3_HTTPS_BASE = "https://czi-biohub-references.s3.us-west-2.amazonaws.com/guidescan-indexes/"

_HEADER = """\
# Reference files configuration for cellxgene-schema-cli
# This file maps reference data keys to their download URLs
# Used by ReferenceFileManager to download and cache reference files
# Auto-updated by update_guidescan_reference_urls.py when gene_assembly.yml changes.
"""


def fasta_url_to_archive_name(fasta_url: str) -> str:
    """Derive the tar.gz archive name from a fasta_url.

    Mirrors the logic in data-pipelines/run_guidescan_index.py:
    Path(filename).stem strips the last extension (.gz).
    """
    filename = fasta_url.rstrip("/").split("/")[-1]
    stem = Path(filename).stem  # strips .gz
    return f"{stem}.tar.gz"


def main() -> None:
    with open(GENE_ASSEMBLY_PATH) as f:
        gene_assembly = yaml.safe_load(f) or {}

    with open(REFERENCE_FILES_PATH) as f:
        reference_files = yaml.safe_load(f) or {}

    guidescan_indexes: dict = reference_files.setdefault("guidescan_indexes", {})
    changed = False

    for species_key, info in gene_assembly.items():
        fasta_url = info.get("fasta_url", "")
        if not fasta_url:
            print(f"WARNING: no fasta_url for '{species_key}', skipping")
            continue

        new_url = S3_HTTPS_BASE + fasta_url_to_archive_name(fasta_url)
        assembly_version = info.get("assembly_version", "")
        description = info.get("description", species_key)
        new_description = f"GuideScan2 index for {description} ({assembly_version})"

        entry = guidescan_indexes.setdefault(species_key, {})

        if entry.get("url") != new_url:
            print(f"  {species_key}: {entry.get('url') or '(none)'} -> {new_url}")
            entry["url"] = new_url
            changed = True

        if entry.get("description") != new_description:
            entry["description"] = new_description
            changed = True

        # Ensure organism_id is present (sourced from gene_assembly.yml).
        organism_id = info.get("organism_id", "")
        if organism_id and entry.get("organism_id") != organism_id:
            entry["organism_id"] = organism_id
            changed = True

    if not changed:
        print("reference_files.yml is already up to date.")
        return

    yaml_str = yaml.dump(reference_files, default_flow_style=False, allow_unicode=True, sort_keys=False)
    with open(REFERENCE_FILES_PATH, "w") as f:
        f.write(_HEADER)
        f.write(yaml_str)

    print("Updated reference_files.yml")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
