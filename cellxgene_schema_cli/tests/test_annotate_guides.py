"""
Tests for guidescan2 gRNA annotation functionality.
"""

import os
from pathlib import Path

import pandas as pd
import pytest
from cellxgene_schema import annotate_guides


# Fixtures
@pytest.fixture
def fixtures_dir():
    """Return path to guidescan2 test fixtures directory."""
    return Path(__file__).parent / "fixtures" / "guidescan2"


@pytest.fixture
def sample_guidescan_csv(fixtures_dir):
    """Return path to sample guidescan2 output CSV."""
    return str(fixtures_dir / "sample_guidescan_output.csv")


@pytest.fixture
def sample_gene_coordinates_csv(fixtures_dir):
    """Return path to sample gene coordinates CSV."""
    return str(fixtures_dir / "sample_gene_coordinates.csv")


@pytest.fixture
def expected_output_csv(fixtures_dir):
    """Return path to expected annotated output CSV."""
    return str(fixtures_dir / "expected_annotated_output.csv")


@pytest.fixture
def sample_gene_coordinates_df(sample_gene_coordinates_csv):
    """Load sample gene coordinates as DataFrame."""
    df = pd.read_csv(sample_gene_coordinates_csv)
    df = df.rename(columns={"chromosome": "chrom"})
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    return df


# Test individual functions


def test_get_best_matches(sample_guidescan_csv):
    """Test that best matches are correctly selected per guide."""
    result = annotate_guides.get_best_matches(sample_guidescan_csv)

    # Should have 6 unique guides
    assert len(result) == 6

    # guide1 should select distance=0 over distance=1
    guide1_rows = result[result["id"] == "guide1"]
    assert len(guide1_rows) == 1
    assert guide1_rows.iloc[0]["distance"] == "0"
    assert guide1_rows.iloc[0]["chromosome"] == "1"

    # guide6 should select distance=0 over distance=2
    guide6_rows = result[result["id"] == "guide6"]
    assert len(guide6_rows) == 1
    assert guide6_rows.iloc[0]["distance"] == "0"
    assert guide6_rows.iloc[0]["position"] == "3005"

    # guide3 should have NA match
    guide3_rows = result[result["id"] == "guide3"]
    assert len(guide3_rows) == 1
    assert guide3_rows.iloc[0]["chromosome"] == "NA"


def test_split_exact_and_no_matches(sample_guidescan_csv):
    """Test splitting of exact matches and no matches."""
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv)
    exact, no_match = annotate_guides.split_exact_and_no_matches(best_matches)

    # Should have 5 exact matches
    assert len(exact) == 5
    assert "guide1" in exact["id"].values
    assert "guide2" in exact["id"].values
    assert "guide4" in exact["id"].values
    assert "guide5" in exact["id"].values
    assert "guide6" in exact["id"].values

    # Should have 1 no match
    assert len(no_match) == 1
    assert "guide3" in no_match["id"].values


def test_format_exact_matches(sample_guidescan_csv):
    """Test formatting of exact matches with coordinate conversion."""
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv)
    exact, _ = annotate_guides.split_exact_and_no_matches(best_matches)
    formatted = annotate_guides.format_exact_matches(exact)

    # Should have same number of rows
    assert len(formatted) == len(exact)

    # Check guide1 (+ strand, position 1000)
    guide1 = formatted[formatted["id"] == "guide1"].iloc[0]
    assert guide1["sequence"] == "TGCCTCGCGCAGCTCGCGG"  # 20bp protospacer
    assert guide1["pam"] == "NGG"  # 3bp PAM
    assert guide1["chromosome"] == "1"
    assert guide1["start"] == 999  # 1000 - 1 (convert to 0-based)
    assert guide1["end"] == 1019  # 1000 + 20 - 1
    assert guide1["sense"] == "+"

    # Check guide2 (- strand, position 1982)
    guide2 = formatted[formatted["id"] == "guide2"].iloc[0]
    assert guide2["sequence"] == "GAGTTCGCTGCGCGCTGTT"
    assert guide2["pam"] == "NGG"
    assert guide2["chromosome"] == "3"
    assert guide2["start"] == 1962  # 1982 - 20
    assert guide2["end"] == 1982  # position
    assert guide2["sense"] == "-"


def test_normalize_chromosomes_guides_with_chr_prefix():
    """Test chromosome normalization when guides have chr prefix but genes don't."""
    guides_df = pd.DataFrame({"id": ["g1"], "chromosome": ["chr1"], "start": [100], "end": [200]})
    genes_df = pd.DataFrame({"chrom": ["1"], "start": [150], "end": [250]})

    guides_norm, genes_norm = annotate_guides.normalize_chromosomes(guides_df, genes_df)

    assert genes_norm["chrom"].iloc[0] == "chr1"


def test_normalize_chromosomes_guides_without_chr_prefix():
    """Test chromosome normalization when guides don't have chr prefix but genes do."""
    guides_df = pd.DataFrame({"id": ["g1"], "chromosome": ["1"], "start": [100], "end": [200]})
    genes_df = pd.DataFrame({"chrom": ["chr1"], "start": [150], "end": [250]})

    guides_norm, genes_norm = annotate_guides.normalize_chromosomes(guides_df, genes_df)

    assert genes_norm["chrom"].iloc[0] == "1"


def test_normalize_chromosomes_both_match():
    """Test chromosome normalization when both use the same format."""
    guides_df = pd.DataFrame({"id": ["g1"], "chromosome": ["1"], "start": [100], "end": [200]})
    genes_df = pd.DataFrame({"chrom": ["1"], "start": [150], "end": [250]})

    guides_norm, genes_norm = annotate_guides.normalize_chromosomes(guides_df, genes_df)

    assert genes_norm["chrom"].iloc[0] == "1"


def test_find_gene_overlaps(sample_guidescan_csv, sample_gene_coordinates_df):
    """Test bioframe overlap detection."""
    # Get formatted guides
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv)
    exact, _ = annotate_guides.split_exact_and_no_matches(best_matches)
    formatted_guides = annotate_guides.format_exact_matches(exact)

    # Find overlaps
    overlaps = annotate_guides.find_gene_overlaps(formatted_guides, sample_gene_coordinates_df)

    # Should find overlaps for guide1 (overlaps 2 genes), guide2, guide5 (overlaps 2 genes), guide6
    # guide4 should have no overlaps (different chromosome region)
    assert len(overlaps) > 0

    # Check that guide1 overlaps both BRCA1 and BRCA1-AS1
    guide1_overlaps = overlaps[overlaps["id_"] == "guide1"]
    assert len(guide1_overlaps) == 2
    gene_names = set(guide1_overlaps["gene_name"].values)
    assert "BRCA1" in gene_names
    assert "BRCA1-AS1" in gene_names


def test_create_annotated_output(sample_guidescan_csv, sample_gene_coordinates_df):
    """Test creation of final annotated output."""
    # Process all steps
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv)
    exact, no_match = annotate_guides.split_exact_and_no_matches(best_matches)
    formatted_guides = annotate_guides.format_exact_matches(exact)
    overlaps = annotate_guides.find_gene_overlaps(formatted_guides, sample_gene_coordinates_df)

    # Create annotated output
    result = annotate_guides.create_annotated_output(overlaps, formatted_guides, no_match)

    # Should have rows for:
    # - guide1 (2 overlaps with genes)
    # - guide2 (1 overlap)
    # - guide3 (1 no-match row)
    # - guide4 (1 no-overlap row)
    # - guide5 (2 overlaps with genes)
    # - guide6 (1 overlap)
    assert len(result) == 8

    # Check guide3 (no genomic match)
    guide3_rows = result[result["id"] == "guide3"]
    assert len(guide3_rows) == 1
    assert guide3_rows.iloc[0]["chromosome"] == "NA"
    assert guide3_rows.iloc[0]["gene_id"] == "NA"
    assert guide3_rows.iloc[0]["gene_name"] == "NA"

    # Check guide4 (genomic match but no gene overlap)
    guide4_rows = result[result["id"] == "guide4"]
    assert len(guide4_rows) == 1
    assert guide4_rows.iloc[0]["chromosome"] == "10"
    assert guide4_rows.iloc[0]["gene_id"] == "NA"
    assert guide4_rows.iloc[0]["gene_name"] == "NA"

    # Check guide1 (overlaps 2 genes)
    guide1_rows = result[result["id"] == "guide1"]
    assert len(guide1_rows) == 2
    gene_names = set(guide1_rows["gene_name"].values)
    assert "BRCA1" in gene_names
    assert "BRCA1-AS1" in gene_names


def test_end_to_end_annotation(sample_guidescan_csv, sample_gene_coordinates_csv, tmp_path):
    """Test full end-to-end annotation pipeline."""
    output_csv = str(tmp_path / "annotated_output.csv")

    # Mock the load_gene_coordinates function to use our test data
    original_load = annotate_guides.load_gene_coordinates

    def mock_load_gene_coordinates(species):
        df = pd.read_csv(sample_gene_coordinates_csv)
        df = df.rename(columns={"chromosome": "chrom"})
        df["start"] = df["start"].astype("int64")
        df["end"] = df["end"].astype("int64")
        return df

    # Temporarily replace the function
    annotate_guides.load_gene_coordinates = mock_load_gene_coordinates

    try:
        # Run annotation
        annotate_guides.annotate_guides_from_guidescan_csv(sample_guidescan_csv, "test_species", output_csv)

        # Check that output file was created
        assert os.path.exists(output_csv)

        # Load and verify output
        output_df = pd.read_csv(output_csv)
        assert len(output_df) == 8

        # Verify required columns exist
        required_cols = ["id", "sequence", "pam", "chromosome", "start", "end", "sense", "gene_id", "gene_name"]
        for col in required_cols:
            assert col in output_df.columns

        # Verify specific guides
        guide1_rows = output_df[output_df["id"] == "guide1"]
        assert len(guide1_rows) == 2  # Overlaps 2 genes

        guide3_rows = output_df[output_df["id"] == "guide3"]
        assert len(guide3_rows) == 1
        assert str(guide3_rows.iloc[0]["chromosome"]) == "NA"

    finally:
        # Restore original function
        annotate_guides.load_gene_coordinates = original_load


def test_coordinate_conversion():
    """Test that coordinate conversion from 1-based to 0-based BED format is correct."""
    # Create test data with known coordinates
    test_df = pd.DataFrame(
        [
            {
                "id": "test1",
                "sequence": "AAAAAAAAAAAAAAAAAAAAAGG",
                "chromosome": "1",
                "position": "1000",
                "strand": "+",
                "distance": "0",
            },
            {
                "id": "test2",
                "sequence": "TTTTTTTTTTTTTTTTTTTTTNGG",
                "chromosome": "1",
                "position": "2000",
                "strand": "-",
                "distance": "0",
            },
        ]
    )

    formatted = annotate_guides.format_exact_matches(test_df)

    # For + strand at position 1000:
    # - guidescan2: 1-based position 1000
    # - BED: start = 999 (0-based), end = 1019 (20bp sequence)
    test1 = formatted[formatted["id"] == "test1"].iloc[0]
    assert test1["start"] == 999
    assert test1["end"] == 1019

    # For - strand at position 2000:
    # - guidescan2: 1-based position 2000 (end of sequence)
    # - BED: start = 1980, end = 2000
    test2 = formatted[formatted["id"] == "test2"].iloc[0]
    assert test2["start"] == 1980
    assert test2["end"] == 2000


def test_handles_empty_input():
    """Test that functions handle empty DataFrames gracefully."""
    empty_df = pd.DataFrame()

    # format_exact_matches should return empty DataFrame with correct columns
    result = annotate_guides.format_exact_matches(empty_df)
    assert len(result) == 0
    assert "id" in result.columns
    assert "sequence" in result.columns

    # find_gene_overlaps should return empty DataFrame
    genes_df = pd.DataFrame({"chrom": [], "start": [], "end": [], "gene_id": [], "gene_name": []})
    result = annotate_guides.find_gene_overlaps(empty_df, genes_df)
    assert len(result) == 0
