#!/usr/bin/env python3
"""Detect which perturbation-support species had their fasta_url changed in gene_info.yml.

Compares the current gene_info.yml against HEAD~1. For each species with
perturbation_support: true, checks whether the fully-expanded fasta_url
(after substituting {version} and {assembly_version}) has changed.

Writes a space-separated list of changed species keys to the GITHUB_OUTPUT
file (if running in GitHub Actions) and also prints them to stdout.

Exit codes:
  0 — ran successfully (no changed species is still a success)
  1 — unexpected error
"""

import os
import subprocess
import sys
from pathlib import Path

import yaml

GENE_INFO_REPO_PATH = "cellxgene_schema_cli/cellxgene_schema/gencode_files/gene_info.yml"


def expand_fasta_url(species_info: dict) -> str:
    """Expand fasta_url template with version and assembly_version."""
    fasta_url = species_info.get("fasta_url", "")
    if not fasta_url:
        return ""
    version = str(species_info.get("version", ""))
    assembly_version = str(species_info.get("assembly_version", ""))
    return fasta_url.format(version=version, assembly_version=assembly_version)


def get_gene_info_at_head_minus_1() -> dict:
    """Return gene_info.yml content from the parent commit (HEAD~1)."""
    result = subprocess.run(
        ["git", "show", f"HEAD~1:{GENE_INFO_REPO_PATH}"],
        capture_output=True,
        text=True,
        check=True,
    )
    return yaml.safe_load(result.stdout) or {}


def get_current_gene_info() -> dict:
    """Return gene_info.yml content from the working tree."""
    repo_root = Path(__file__).parents[2]
    path = repo_root / GENE_INFO_REPO_PATH
    with open(path) as f:
        return yaml.safe_load(f) or {}


def detect_changed_species() -> list[str]:
    """Return species keys whose expanded fasta_url changed since HEAD~1."""
    try:
        prev = get_gene_info_at_head_minus_1()
    except subprocess.CalledProcessError:
        # No parent commit (e.g. initial commit) — treat everything as new.
        prev = {}

    curr = get_current_gene_info()

    changed = []
    for key, info in curr.items():
        if not info.get("perturbation_support"):
            continue
        if not info.get("fasta_url"):
            continue

        curr_url = expand_fasta_url(info)
        prev_url = expand_fasta_url(prev.get(key, {}))

        if curr_url != prev_url:
            changed.append(key)
            print(f"  {key}: fasta_url changed")
            if prev_url:
                print(f"    was: {prev_url}")
            print(f"    now: {curr_url}")

    return changed


def main() -> None:
    print("Detecting fasta_url changes for perturbation-support species...")
    changed = detect_changed_species()

    changed_str = " ".join(changed)
    if changed:
        print(f"\nChanged species: {changed_str}")
    else:
        print("\nNo fasta_url changes detected.")

    github_output = os.environ.get("GITHUB_OUTPUT")
    if github_output:
        with open(github_output, "a") as f:
            f.write(f"changed_species={changed_str}\n")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
