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


def test_update_h5ad_with_guide_annotations():
    """Test updating h5ad with guide annotations."""
    import anndata as ad

    # Create a minimal AnnData object with genetic perturbations
    adata = ad.AnnData()
    adata.uns["genetic_perturbations"] = {
        "guide1": {
            "target_genomic_regions": [],
            "target_features": {},
        },
        "guide2": {
            "target_genomic_regions": [],
            "target_features": {},
        },
    }

    # Create annotations DataFrame with multiple gene overlaps
    annotations_df = pd.DataFrame(
        [
            {
                "id": "guide1",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG00000123456",
                "gene_name": "BRCA1",
            },
            {
                "id": "guide1",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG00000234567",
                "gene_name": "BRCA1-AS1",
            },
            {
                "id": "guide2",
                "chromosome": "3",
                "start": 1963,
                "end": 1982,
                "sense": "-",
                "gene_id": "ENSG00000345678",
                "gene_name": "EGFR",
            },
        ]
    )

    # Update the h5ad
    result_adata = annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)

    # Check guide1 (overlaps 2 genes at same location)
    guide1 = result_adata.uns["genetic_perturbations"]["guide1"]
    # Only genomic annotations should be updated (not sequence/PAM)
    assert len(guide1["target_genomic_regions"]) == 1
    # Schema 7.1.0: 1-based coordinates (BED 999-1019 becomes 1000-1019)
    assert "1:1000-1019(+)" in guide1["target_genomic_regions"]
    assert len(guide1["target_features"]) == 2
    assert guide1["target_features"]["ENSG00000123456"] == "BRCA1"
    assert guide1["target_features"]["ENSG00000234567"] == "BRCA1-AS1"

    # Check guide2 (overlaps 1 gene)
    guide2 = result_adata.uns["genetic_perturbations"]["guide2"]
    # Only genomic annotations should be updated (not sequence/PAM)
    assert len(guide2["target_genomic_regions"]) == 1
    # Schema 7.1.0: 1-based coordinates (BED 1963-1982 becomes 1964-1982)
    assert "3:1964-1982(-)" in guide2["target_genomic_regions"]
    assert len(guide2["target_features"]) == 1
    assert guide2["target_features"]["ENSG00000345678"] == "EGFR"


def test_update_h5ad_with_guide_annotations_skips_na():
    """Test that guides with no gene overlaps are skipped."""
    import anndata as ad

    adata = ad.AnnData()
    adata.uns["genetic_perturbations"] = {
        "guide1": {
            "target_genomic_regions": ["old:1-10(+)"],
            "target_features": {"OLD_GENE": "OLD_NAME"},
        }
    }

    # Annotations with no gene overlaps (all NA)
    annotations_df = pd.DataFrame(
        [
            {
                "id": "guide1",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "NA",
                "gene_name": "NA",
            }
        ]
    )

    # Update should skip this guide
    result_adata = annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)

    # All values should be unchanged (guide was skipped due to no gene overlaps)
    guide1 = result_adata.uns["genetic_perturbations"]["guide1"]
    assert guide1["target_genomic_regions"] == ["old:1-10(+)"]
    assert guide1["target_features"] == {"OLD_GENE": "OLD_NAME"}


def test_update_h5ad_with_guide_annotations_missing_guide_error():
    """Test that missing guide IDs raise ValueError."""
    import anndata as ad

    adata = ad.AnnData()
    adata.uns["genetic_perturbations"] = {
        "guide1": {},
    }

    # Annotations for guide that doesn't exist
    annotations_df = pd.DataFrame(
        [
            {
                "id": "guide_missing",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG123",
                "gene_name": "GENE1",
            }
        ]
    )

    # Should raise ValueError
    with pytest.raises(ValueError, match="missing from adata.uns"):
        annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)


def test_update_h5ad_with_guide_annotations_no_genetic_perturbations():
    """Test that missing genetic_perturbations raises KeyError."""
    import anndata as ad

    adata = ad.AnnData()
    # No genetic_perturbations in uns

    annotations_df = pd.DataFrame(
        [
            {
                "id": "guide1",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG123",
                "gene_name": "GENE1",
            }
        ]
    )

    # Should raise KeyError
    with pytest.raises(KeyError, match="genetic_perturbations"):
        annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)


def test_update_h5ad_chromosome_ensembl_format():
    """Test that chromosomes are converted to ENSEMBL format per schema 7.1.0."""
    import anndata as ad

    adata = ad.AnnData()
    adata.uns["genetic_perturbations"] = {
        "guide1": {},
        "guide2": {},
        "guide3": {},
    }

    # Test various chromosome formats
    annotations_df = pd.DataFrame(
        [
            # chr prefix should be removed
            {
                "id": "guide1",
                "chromosome": "chr1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG123",
                "gene_name": "GENE1",
            },
            # chrM should become MT
            {
                "id": "guide2",
                "chromosome": "chrM",
                "start": 100,
                "end": 120,
                "sense": "+",
                "gene_id": "ENSG456",
                "gene_name": "GENE2",
            },
            # M should become MT
            {
                "id": "guide3",
                "chromosome": "M",
                "start": 200,
                "end": 220,
                "sense": "-",
                "gene_id": "ENSG789",
                "gene_name": "GENE3",
            },
        ]
    )

    result_adata = annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)

    # chr1 should be converted to 1
    guide1 = result_adata.uns["genetic_perturbations"]["guide1"]
    assert "1:1000-1019(+)" in guide1["target_genomic_regions"]

    # chrM should be converted to MT
    guide2 = result_adata.uns["genetic_perturbations"]["guide2"]
    assert "MT:101-120(+)" in guide2["target_genomic_regions"]

    # M should be converted to MT
    guide3 = result_adata.uns["genetic_perturbations"]["guide3"]
    assert "MT:201-220(-)" in guide3["target_genomic_regions"]


def test_update_h5ad_removes_ensembl_version():
    """Test that Ensembl ID version numbers are removed per schema 7.1.0."""
    import anndata as ad

    adata = ad.AnnData()
    adata.uns["genetic_perturbations"] = {
        "guide1": {},
        "guide2": {},
    }

    # Test gene IDs with and without versions
    annotations_df = pd.DataFrame(
        [
            # Ensembl ID with version number
            {
                "id": "guide1",
                "chromosome": "1",
                "start": 999,
                "end": 1019,
                "sense": "+",
                "gene_id": "ENSG00000186092.7",
                "gene_name": "GENE1",
            },
            # Non-Ensembl ID should not be modified
            {
                "id": "guide2",
                "chromosome": "2",
                "start": 100,
                "end": 120,
                "sense": "+",
                "gene_id": "SOME_OTHER_ID.5",
                "gene_name": "GENE2",
            },
        ]
    )

    result_adata = annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)

    # Ensembl ID version should be removed
    guide1 = result_adata.uns["genetic_perturbations"]["guide1"]
    assert "ENSG00000186092" in guide1["target_features"]
    assert "ENSG00000186092.7" not in guide1["target_features"]
    assert guide1["target_features"]["ENSG00000186092"] == "GENE1"

    # Non-Ensembl ID should keep version
    guide2 = result_adata.uns["genetic_perturbations"]["guide2"]
    assert "SOME_OTHER_ID.5" in guide2["target_features"]
    assert guide2["target_features"]["SOME_OTHER_ID.5"] == "GENE2"
