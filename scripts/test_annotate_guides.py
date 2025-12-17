#!/usr/bin/env python3
"""Test script to validate annotate_guides.py with real dataset.

This script:
1. Loads small_perturbations_cleaned.h5ad (input with guide sequences but no annotations)
2. Runs annotate_perturbations_in_h5ad() to add genomic annotations
3. Compares the result to small_perturbations.h5ad (expected output)
4. Reports differences and validates the annotation process
"""

import sys
from pathlib import Path

import anndata as ad

# Add parent directory to path to import cellxgene_schema
script_dir = Path(__file__).parent.absolute()
repo_root = script_dir.parent
sys.path.insert(0, str(repo_root / "cellxgene_schema_cli"))

from cellxgene_schema.annotate_guides import annotate_perturbations_in_h5ad


def compare_genetic_perturbations(adata1: ad.AnnData, adata2: ad.AnnData, label1: str, label2: str) -> None:
    """Compare genetic_perturbations between two AnnData objects.

    :param adata1: First AnnData object
    :param adata2: Second AnnData object
    :param label1: Label for first dataset (for display)
    :param label2: Label for second dataset (for display)
    """
    print(f"\n{'='*80}")
    print(f"Comparing genetic_perturbations: {label1} vs {label2}")
    print(f"{'='*80}\n")

    # Check if both have genetic_perturbations
    has_gp1 = "genetic_perturbations" in adata1.uns
    has_gp2 = "genetic_perturbations" in adata2.uns

    if not has_gp1 and not has_gp2:
        print("✓ Neither dataset has genetic_perturbations")
        return
    elif not has_gp1:
        print(f"✗ {label1} missing genetic_perturbations")
        return
    elif not has_gp2:
        print(f"✗ {label2} missing genetic_perturbations")
        return

    gp1 = adata1.uns["genetic_perturbations"]
    gp2 = adata2.uns["genetic_perturbations"]

    # Compare guide IDs
    guides1 = set(gp1.keys())
    guides2 = set(gp2.keys())

    print("Guide counts:")
    print(f"  {label1}: {len(guides1)} guides")
    print(f"  {label2}: {len(guides2)} guides")

    if guides1 != guides2:
        print("\n✗ Guide ID mismatch:")
        only_in_1 = guides1 - guides2
        only_in_2 = guides2 - guides1
        if only_in_1:
            print(f"  Only in {label1}: {sorted(only_in_1)}")
        if only_in_2:
            print(f"  Only in {label2}: {sorted(only_in_2)}")
    else:
        print(f"\n✓ Both datasets have the same {len(guides1)} guide IDs")

    # Compare each guide's annotations
    print(f"\n{'─'*80}")
    print("Per-guide comparison:")
    print(f"{'─'*80}\n")

    common_guides = sorted(guides1 & guides2)
    differences_found = False

    for guide_id in common_guides:
        guide1 = gp1[guide_id]
        guide2 = gp2[guide_id]

        guide_diffs = []

        # Compare all fields (skip target_features since it's not in expected output)
        fields = set(guide1.keys()) | set(guide2.keys())
        fields_to_compare = fields - {"target_features"}  # Skip target_features for now

        for field in sorted(fields_to_compare):
            val1 = guide1.get(field)
            val2 = guide2.get(field)

            # Compare values (handle lists/dicts specially)
            if isinstance(val1, (list, dict)) or isinstance(val2, (list, dict)):
                if str(val1) != str(val2):
                    guide_diffs.append((field, val1, val2))
            else:
                if val1 != val2:
                    guide_diffs.append((field, val1, val2))

        # Report differences for this guide
        if guide_diffs:
            differences_found = True
            print(f"\nGuide: {guide_id}")
            for field, val1, val2 in guide_diffs:
                print(f"  Field: {field}")
                print(f"    {label1}: {val1}")
                print(f"    {label2}: {val2}")
        else:
            print(f"✓ Guide {guide_id}: identical")

    if not differences_found:
        print(f"\n{'='*80}")
        print("✓ All guides have identical annotations!")
        print(f"{'='*80}")
    else:
        print(f"\n{'='*80}")
        print("✗ Differences found in guide annotations")
        print(f"{'='*80}")


def main() -> None:
    """Main test script."""
    # File paths
    input_file = script_dir / "small_perturbations_cleaned.h5ad"
    expected_file = script_dir / "small_perturbations.h5ad"

    print("=" * 80)
    print("Testing annotate_guides.py with real dataset")
    print("=" * 80)

    # Check files exist
    if not input_file.exists():
        print(f"✗ Input file not found: {input_file}")
        sys.exit(1)

    if not expected_file.exists():
        print(f"✗ Expected file not found: {expected_file}")
        sys.exit(1)

    print(f"\nInput file:    {input_file}")
    print(f"Expected file: {expected_file}")

    # Load input file
    print("\n" + "─" * 80)
    print("Loading input file...")
    print("─" * 80)
    adata_input = ad.read_h5ad(input_file)
    print(f"✓ Loaded {input_file.name}")
    print(f"  Shape: {adata_input.shape}")
    print(f"  Organism: {adata_input.uns.get('organism_ontology_term_id', 'N/A')}")

    # Check if genetic_perturbations exists
    if "genetic_perturbations" not in adata_input.uns:
        print("\n✗ Input file has no genetic_perturbations in uns!")
        sys.exit(1)

    print(f"  Genetic perturbations: {len(adata_input.uns['genetic_perturbations'])} guides")

    # Show sample guide before annotation
    sample_guide_id = list(adata_input.uns["genetic_perturbations"].keys())[0]
    sample_guide = adata_input.uns["genetic_perturbations"][sample_guide_id]
    print(f"\nSample guide before annotation ({sample_guide_id}):")
    for key, val in sample_guide.items():
        print(f"  {key}: {val}")

    # Run annotation
    print("\n" + "─" * 80)
    print("Running annotation pipeline...")
    print("─" * 80)

    try:
        adata_annotated = annotate_perturbations_in_h5ad(adata_input)
        print("✓ Annotation completed successfully")
        print("\nChecking annotation results...")
        # Count how many guides have target_genomic_regions
        guides_with_regions = 0
        guides_with_features = 0
        for _, guide_data in adata_annotated.uns["genetic_perturbations"].items():
            if "target_genomic_regions" in guide_data and guide_data["target_genomic_regions"]:
                guides_with_regions += 1
            if "target_features" in guide_data and guide_data["target_features"]:
                guides_with_features += 1
        print(
            f"  Guides with target_genomic_regions: {guides_with_regions}/{len(adata_annotated.uns['genetic_perturbations'])}"
        )
        print(
            f"  Guides with target_features: {guides_with_features}/{len(adata_annotated.uns['genetic_perturbations'])}"
        )
    except Exception as e:
        print(f"✗ Annotation failed: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)

    # Show sample guide after annotation
    sample_guide_annotated = adata_annotated.uns["genetic_perturbations"][sample_guide_id]
    print(f"\nSample guide after annotation ({sample_guide_id}):")
    for key, val in sample_guide_annotated.items():
        print(f"  {key}: {val}")

    # Load expected output
    print("\n" + "─" * 80)
    print("Loading expected output file...")
    print("─" * 80)
    adata_expected = ad.read_h5ad(expected_file)
    print(f"✓ Loaded {expected_file.name}")
    print(f"  Shape: {adata_expected.shape}")

    if "genetic_perturbations" in adata_expected.uns:
        print(f"  Genetic perturbations: {len(adata_expected.uns['genetic_perturbations'])} guides")

    # Compare results
    compare_genetic_perturbations(adata_annotated, adata_expected, "Annotated", "Expected")

    print("\n" + "=" * 80)
    print("Test complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
