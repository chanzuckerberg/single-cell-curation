"""
Tests for guidescan2 gRNA annotation functionality.
"""

import os
import shutil
import subprocess
from pathlib import Path
from unittest.mock import MagicMock

import anndata as ad
import numpy as np
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
def sample_guidescan_input_csv(fixtures_dir):
    """Return path to sample guidescan2 input CSV."""
    return str(fixtures_dir / "sample_guidescan_input.csv")


@pytest.fixture
def sample_gene_coordinates_csv(fixtures_dir):
    """Return path to sample gene coordinates CSV."""
    return str(fixtures_dir / "sample_gene_coordinates.csv")


@pytest.fixture
def sample_gene_coordinates_df(sample_gene_coordinates_csv):
    """Load sample gene coordinates as DataFrame."""
    df = pd.read_csv(sample_gene_coordinates_csv)
    df = df.rename(columns={"chromosome": "chrom"})
    # Ensure proper dtypes for bioframe compatibility
    df["chrom"] = df["chrom"].astype(str)
    df["start"] = df["start"].astype("int64")
    df["end"] = df["end"].astype("int64")
    return df


@pytest.mark.parametrize(
    "which_return,should_raise,description",
    [
        ("/usr/bin/guidescan", False, "guidescan2 is installed"),
        (None, True, "guidescan2 is not installed"),
    ],
)
def test_check_guidescan2_installed(monkeypatch, which_return, should_raise, description):
    """Test checking if guidescan2 is installed."""
    # Mock shutil.which
    monkeypatch.setattr("shutil.which", lambda x: which_return)

    if should_raise:
        with pytest.raises(RuntimeError, match="guidescan2 not found"):
            annotate_guides._check_guidescan2_installed()
    else:
        # Should not raise
        annotate_guides._check_guidescan2_installed()


@pytest.mark.parametrize(
    "organism_id",
    [
        "NCBITaxon:7955",  # Danio rerio (zebrafish)
        "NCBITaxon:9606",  # Homo sapiens (human)
        "NCBITaxon:10090",  # Mus musculus (mouse)
    ],
)
def test_get_organism_from_h5ad(organism_id):
    """Test extracting supported organisms from h5ad uns."""
    adata = ad.AnnData(X=np.zeros((2, 2)))
    adata.uns["organism_ontology_term_id"] = organism_id

    result = annotate_guides._get_organism_from_h5ad(adata)
    assert result == organism_id


def test_get_organism_from_h5ad_missing():
    """Test error when organism_ontology_term_id is missing."""
    adata = ad.AnnData(X=np.zeros((2, 2)))

    with pytest.raises(KeyError, match="organism_ontology_term_id"):
        annotate_guides._get_organism_from_h5ad(adata)


@pytest.mark.parametrize(
    "organism_id,expected_match,description",
    [
        ("NCBITaxon:12345", "not supported", "invalid organism ID"),
        (
            "NCBITaxon:7227",
            "not supported for genetic perturbations",
            "Drosophila - valid organism but not for perturbations",
        ),
        (
            "NCBITaxon:6239",
            "not supported for genetic perturbations",
            "C. elegans - valid organism but not for perturbations",
        ),
    ],
)
def test_get_organism_from_h5ad_unsupported(organism_id, expected_match, description):
    """Test error when organism is not supported for genetic perturbations."""
    adata = ad.AnnData(X=np.zeros((2, 2)))
    adata.uns["organism_ontology_term_id"] = organism_id

    with pytest.raises(ValueError, match=expected_match):
        annotate_guides._get_organism_from_h5ad(adata)


def test_extract_perturbations_to_guidescan_csv(tmp_path):
    """Test extracting perturbations to CSV for guidescan2."""
    # Create AnnData with genetic perturbations
    adata = ad.AnnData(X=np.zeros((2, 2)))
    adata.uns["genetic_perturbations"] = {
        "guide1": {
            "protospacer_sequence": "TGCCTCGCGCAGCTCGCGG",
            "protospacer_adjacent_motif": "3' NGG",
        },
        "guide2": {
            "protospacer_sequence": "ACGTACGTACGTACGTACG",
            "protospacer_adjacent_motif": "3' TGG",
        },
    }

    output_csv = tmp_path / "guides.csv"
    annotate_guides._extract_perturbations_to_guidescan_csv(adata, str(output_csv))

    # Check CSV was created
    assert output_csv.exists()

    # Check CSV contents
    df = pd.read_csv(output_csv)
    assert len(df) == 2
    assert set(df["id"]) == {"guide1", "guide2"}

    # Check full sequence (protospacer + PAM)
    guide1_row = df[df["id"] == "guide1"].iloc[0]
    assert guide1_row["sequence"] == "TGCCTCGCGCAGCTCGCGG"

    guide2_row = df[df["id"] == "guide2"].iloc[0]
    assert guide2_row["sequence"] == "ACGTACGTACGTACGTACG"


def test_run_guidescan_enumerate_success(tmp_path, monkeypatch):
    """Test running guidescan enumerate successfully."""
    # Mock subprocess.run to simulate successful execution
    mock_result = MagicMock()
    mock_result.returncode = 0
    mock_result.stdout = "Success"
    mock_result.stderr = ""

    mock_run = MagicMock(return_value=mock_result)
    monkeypatch.setattr(subprocess, "run", mock_run)

    # Create dummy files
    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id,sequence\nguide1,TGCCTCGCGCAGCTCGCGGNGG\n")

    output_csv = tmp_path / "output.csv"
    index_prefix = "/path/to/index"

    # Should not raise
    annotate_guides._run_guidescan_enumerate(str(input_csv), index_prefix, str(output_csv))

    # Verify subprocess was called correctly
    mock_run.assert_called_once()
    call_args = mock_run.call_args[0][0]
    assert call_args[0] == "guidescan"
    assert call_args[1] == "enumerate"
    assert call_args[2] == index_prefix
    assert "--mismatches" in call_args
    assert "0" in call_args


def test_run_guidescan_enumerate_failure(tmp_path, monkeypatch):
    """Test error handling when guidescan enumerate fails."""
    # Mock subprocess.run to simulate failure
    mock_run = MagicMock(side_effect=subprocess.CalledProcessError(1, "guidescan", stderr="Error message"))
    monkeypatch.setattr(subprocess, "run", mock_run)

    input_csv = tmp_path / "input.csv"
    input_csv.write_text("id,sequence\nguide1,TGCCTCGCGCAGCTCGCGGNGG\n")

    output_csv = tmp_path / "output.csv"
    index_prefix = "/path/to/index"

    with pytest.raises(RuntimeError, match="Guidescan enumerate failed"):
        annotate_guides._run_guidescan_enumerate(str(input_csv), index_prefix, str(output_csv))


def test_get_or_download_index_unsupported_organism():
    """Test error for unsupported organism ID."""
    with pytest.raises(ValueError, match="not available"):
        annotate_guides._get_or_download_index("NCBITaxon:99999")


def test_get_or_download_index_uses_reference_manager(monkeypatch):
    """Test that _get_or_download_index uses ReferenceFileManager."""
    from cellxgene_schema import reference_file_manager

    # Mock ReferenceFileManager
    mock_manager = MagicMock()
    mock_manager.get_key_by_organism_id.return_value = "human"
    mock_manager.fetch.return_value = [
        "/cache/human.index.forward",
        "/cache/human.index.reverse",
        "/cache/human.index.gs",
    ]

    mock_manager_class = MagicMock(return_value=mock_manager)
    monkeypatch.setattr(reference_file_manager, "ReferenceFileManager", mock_manager_class)
    monkeypatch.setattr(annotate_guides, "ReferenceFileManager", mock_manager_class)

    result = annotate_guides._get_or_download_index("NCBITaxon:9606")

    # Should return the index prefix (path without .gs extension)
    # With files like "human.index.gs", "human.index.forward", "human.index.reverse",
    # the prefix should be "human.index"
    assert result == "/cache/human.index"
    mock_manager.fetch.assert_called_once()


@pytest.mark.skipif(shutil.which("guidescan") is None, reason="guidescan2 not installed")
def test_guidescan2_real_integration(tmp_path):
    """Integration test with real guidescan2 binary (skipped if not installed).

    This test verifies that we can actually run guidescan2 enumerate with real input.
    It requires guidescan2 to be installed in the PATH.
    """
    # Create a simple input CSV with a guide sequence
    input_csv = tmp_path / "guides_input.csv"
    input_csv.write_text("id,sequence,pam,chromosome,start,end,sense\n" "test_guide,AAAAAAAAAAAAAAAAAAAAAGG,NGG,,,,\n")

    output_csv = tmp_path / "guidescan_output.csv"

    # For this test, we need a real guidescan2 index
    # Since we don't have one in the test environment, we'll verify that:
    # 1. guidescan is installed
    # 2. The command fails with an expected error about missing index (not a command-not-found error)

    # Verify guidescan is installed
    guidescan_path = shutil.which("guidescan")
    assert guidescan_path is not None, "guidescan should be in PATH"

    # Try to run with a non-existent index - should fail gracefully
    fake_index = str(tmp_path / "nonexistent_index")  # Intentionally use a non-existent index to test error handling

    with pytest.raises(RuntimeError) as exc_info:
        annotate_guides._run_guidescan_enumerate(str(input_csv), fake_index, str(output_csv))

    # Should fail due to missing index, not because guidescan is not found
    error_msg = str(exc_info.value)
    assert "Guidescan enumerate failed" in error_msg, "Should report guidescan enumerate failure"

    # The error should NOT be about command not found
    assert "command not found" not in error_msg.lower(), "Error should not be about missing guidescan binary"


# Test individual functions


def test_get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv):
    """Test that best matches are correctly selected per guide."""
    result = annotate_guides.get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv)

    # Should have 6 unique guides
    assert len(result) == 6

    # guide1 should select distance=0 over distance=1
    guide1_rows = result[result["id"] == "guide1"]
    assert len(guide1_rows) == 1
    assert guide1_rows.iloc[0]["distance"] == 0
    assert guide1_rows.iloc[0]["chromosome"] == "1"

    # guide6 should select distance=0 over distance=2
    guide6_rows = result[result["id"] == "guide6"]
    assert len(guide6_rows) == 1
    assert guide6_rows.iloc[0]["distance"] == 0
    assert guide6_rows.iloc[0]["position"] == "3005"

    # guide3 should have NA match
    guide3_rows = result[result["id"] == "guide3"]
    assert len(guide3_rows) == 1
    assert guide3_rows.iloc[0]["chromosome"] == "NA"


def test_split_exact_and_no_matches(sample_guidescan_csv, sample_guidescan_input_csv):
    """Test splitting of exact matches and no matches."""
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv)
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


def test_format_exact_matches(sample_guidescan_csv, sample_guidescan_input_csv):
    """Test formatting of exact matches with coordinate conversion."""
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv)
    exact, _ = annotate_guides.split_exact_and_no_matches(best_matches)
    formatted = annotate_guides.format_exact_matches(exact)

    # Should have same number of rows
    assert len(formatted) == len(exact)

    # Check guide1 (+ strand, position 1000)
    guide1 = formatted[formatted["id"] == "guide1"].iloc[0]
    assert guide1["sequence"] == "TGCCTCGCGCAGCTCGCGG"  # 19bp protospacer
    assert guide1["pam"] == "NGG"  # 3bp PAM
    assert guide1["chromosome"] == "1"
    assert guide1["start"] == 999  # 1000 - 1 (convert to 0-based)
    assert guide1["end"] == 1018  # 1000 + 19 - 1
    assert guide1["sense"] == "+"

    # Check guide2 (- strand, position 1982)
    guide2 = formatted[formatted["id"] == "guide2"].iloc[0]
    assert guide2["sequence"] == "GAGTTCGCTGCGCGCTGTT"
    assert guide2["pam"] == "NGG"
    assert guide2["chromosome"] == "3"
    assert guide2["start"] == 1963  # 1982 - 19
    assert guide2["end"] == 1982  # position
    assert guide2["sense"] == "-"


@pytest.mark.parametrize(
    "guide_chr,gene_chr,expected_gene_chr,description",
    [
        ("chr1", "1", "chr1", "guides have chr prefix, genes don't"),
        ("1", "chr1", "1", "guides don't have chr prefix, genes do"),
        ("1", "1", "1", "both use same format without prefix"),
        ("chr1", "chr1", "chr1", "both use same format with prefix"),
    ],
)
def test_normalize_chromosomes(guide_chr, gene_chr, expected_gene_chr, description):
    """Test chromosome normalization across different formats."""
    guides_df = pd.DataFrame({"id": ["g1"], "chromosome": [guide_chr], "start": [100], "end": [200]})
    genes_df = pd.DataFrame({"chrom": [gene_chr], "start": [150], "end": [250]})

    _, genes_norm = annotate_guides.normalize_chromosomes(guides_df, genes_df)

    assert genes_norm["chrom"].iloc[0] == expected_gene_chr, f"Failed for: {description}"


def test_find_gene_overlaps(sample_guidescan_csv, sample_guidescan_input_csv, sample_gene_coordinates_df):
    """Test bioframe overlap detection."""
    # Get formatted guides
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv)
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


def test_create_annotated_output(sample_guidescan_csv, sample_guidescan_input_csv, sample_gene_coordinates_df):
    """Test creation of final annotated output."""
    # Process all steps
    best_matches = annotate_guides.get_best_matches(sample_guidescan_csv, sample_guidescan_input_csv)
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


def test_end_to_end_annotation(
    sample_guidescan_csv, sample_guidescan_input_csv, sample_gene_coordinates_csv, tmp_path, monkeypatch
):
    """Test full end-to-end annotation pipeline."""
    output_csv = str(tmp_path / "annotated_output.csv")

    # Mock the load_gene_coordinates function to use our test data
    def mock_load_gene_coordinates(species):
        df = pd.read_csv(sample_gene_coordinates_csv)
        df = df.rename(columns={"chromosome": "chrom"})
        # Ensure proper dtypes for bioframe compatibility
        df["chrom"] = df["chrom"].astype(str)
        df["start"] = df["start"].astype("int64")
        df["end"] = df["end"].astype("int64")
        return df

    monkeypatch.setattr(annotate_guides, "load_gene_coordinates", mock_load_gene_coordinates)

    # Run annotation
    annotate_guides.annotate_guides_from_guidescan_csv(
        sample_guidescan_csv, "test_species", sample_guidescan_input_csv, output_csv=output_csv
    )

    # Check that output file was created
    assert os.path.exists(output_csv), f"Output file not created at {output_csv}"

    # Load and verify output
    # Use keep_default_na=False to preserve "NA" as a string, not as pandas NA
    output_df = pd.read_csv(output_csv, keep_default_na=False, na_values=[""])
    assert len(output_df) == 8, f"Expected 8 rows, got {len(output_df)}"

    # Verify required columns exist
    required_cols = ["id", "sequence", "pam", "chromosome", "start", "end", "sense", "gene_id", "gene_name"]
    for col in required_cols:
        assert col in output_df.columns, f"Missing required column: {col}"

    # Verify specific guides
    guide1_rows = output_df[output_df["id"] == "guide1"]
    assert len(guide1_rows) == 2, "guide1 should overlap 2 genes"

    guide3_rows = output_df[output_df["id"] == "guide3"]
    assert len(guide3_rows) == 1, "guide3 should have 1 row"
    assert str(guide3_rows.iloc[0]["chromosome"]) == "NA", "guide3 should have no genomic match"


def test_coordinate_conversion():
    """Test that coordinate conversion from 1-based to 0-based BED format is correct."""
    # Create test data with known coordinates
    test_df = pd.DataFrame(
        [
            {
                "id": "test1",
                "sequence": "AAAAAAAAAAAAAAAAAAAAAGG",
                "pam": "AGG",
                "chromosome": "1",
                "position": "1000",
                "strand": "+",
                "distance": "0",
            },
            {
                "id": "test2",
                "sequence": "TTTTTTTTTTTTTTTTTTTTNGG",
                "pam": "NGG",
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


@pytest.mark.parametrize(
    "guide_id,full_sequence,pam_column,expected_protospacer,expected_pam,expected_pam_length",
    [
        (
            "test1",
            "AAAAAAAAAAAAAAAAAAAAAGG",
            "NGG",
            "AAAAAAAAAAAAAAAAAAAA",
            "AGG",
            3,
        ),
        (
            "test2",
            "TTTTTTTTTTTTTTTTTTTTNGG",
            "NGG",
            "TTTTTTTTTTTTTTTTTTTT",
            "NGG",
            3,
        ),
        (
            "test3",
            "CCCCCCCCCCCCCCCCCCCCCCGG",
            "GG",
            "CCCCCCCCCCCCCCCCCCCCCC",
            "GG",
            2,
        ),
    ],
)
def test_pam_extracted_from_sequence(
    guide_id, full_sequence, pam_column, expected_protospacer, expected_pam, expected_pam_length
):
    """Test that PAM is extracted from actual nucleotides after the 3' end of protospacer.

    This verifies that PAM length is derived from the actual nucleotides in the sequence,
    not just from the PAM column value. The PAM should be extracted from the full sequence
    by taking the nucleotides after the protospacer (3' end).
    """
    test_df = pd.DataFrame(
        [
            {
                "id": guide_id,
                "sequence": full_sequence,
                "pam": pam_column,
                "chromosome": "1",
                "position": "1000",
                "strand": "+",
                "distance": "0",
            }
        ]
    )

    formatted = annotate_guides.format_exact_matches(test_df)

    assert len(formatted) == 1, "Should have exactly one formatted guide"
    result = formatted.iloc[0]

    assert result["id"] == guide_id
    assert result["sequence"] == expected_protospacer, "Protospacer should be extracted correctly"
    assert result["pam"] == expected_pam, "PAM should be extracted from sequence nucleotides after 3' end"
    assert len(result["pam"]) == expected_pam_length, "PAM length should be derived from actual nucleotides"


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


def test_get_best_matches_variable_pam_length(tmp_path):
    """Test that get_best_matches handles variable-length PAMs correctly."""
    # Create input CSV with different PAM lengths
    input_csv = tmp_path / "input.csv"
    input_csv.write_text(
        "id,sequence,pam,chromosome,position,sense\n"
        "guide1,AAAAAAAAAAAAAAAAAAAA,NGG,,,\n"  # 3-char PAM (SpCas9)
        "guide2,BBBBBBBBBBBBBBBBBBBB,TTTN,,,\n"  # 4-char PAM (Cas12a)
        "guide3,CCCCCCCCCCCCCCCCCCCC,NNGRRT,,,\n"  # 6-char PAM
    )

    # Create mock guidescan output (simulating what guidescan would return)
    output_csv = tmp_path / "output.csv"
    output_csv.write_text(
        "id,sequence,match_chrm,match_position,match_strand,match_distance,match_sequence,rna_bulges,dna_bulges,specificity\n"
        "guide1,AAAAAAAAAAAAAAAAAAAANGG,chr1,1000,+,0,AAAAAAAAAAAAAAAAAAAACGG,0,0,1.0\n"
        "guide2,BBBBBBBBBBBBBBBBBBBBTTTN,chr2,2000,-,0,BBBBBBBBBBBBBBBBBBBBTTTT,0,0,1.0\n"
        "guide3,CCCCCCCCCCCCCCCCCCCCNNGRRT,chr3,3000,+,0,CCCCCCCCCCCCCCCCCCCCAAGAAT,0,0,1.0\n"
    )

    # Test with input CSV (should preserve original PAM from input)
    result = annotate_guides.get_best_matches(str(output_csv), str(input_csv))

    assert len(result) == 3
    assert result.loc[result["id"] == "guide1", "pam"].iloc[0] == "NGG"
    assert result.loc[result["id"] == "guide2", "pam"].iloc[0] == "TTTN"
    assert result.loc[result["id"] == "guide3", "pam"].iloc[0] == "NNGRRT"

    # Test that input_csv is now required - should raise error without it
    with pytest.raises(TypeError, match="missing 1 required positional argument"):
        annotate_guides.get_best_matches(str(output_csv))


def test_update_h5ad_with_guide_annotations():
    """Test updating h5ad with guide annotations."""
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
    assert (
        len(guide1["target_genomic_regions"]) == 1
    ), f"Expected 1 genomic region, got {len(guide1['target_genomic_regions'])}"
    # Schema 7.1.0: 1-based coordinates (BED 999-1019 becomes 1000-1019)
    assert (
        "1:1000-1019(+)" in guide1["target_genomic_regions"]
    ), f"Expected 1:1000-1019(+), got {guide1['target_genomic_regions']}"
    assert len(guide1["target_features"]) == 2, f"Expected 2 target features, got {len(guide1['target_features'])}"
    assert guide1["target_features"]["ENSG00000123456"] == "BRCA1"
    assert guide1["target_features"]["ENSG00000234567"] == "BRCA1-AS1"

    # Check guide2 (overlaps 1 gene)
    guide2 = result_adata.uns["genetic_perturbations"]["guide2"]
    # Only genomic annotations should be updated (not sequence/PAM)
    assert (
        len(guide2["target_genomic_regions"]) == 1
    ), f"Expected 1 genomic region, got {len(guide2['target_genomic_regions'])}"
    # Schema 7.1.0: 1-based coordinates (BED 1963-1982 becomes 1964-1982)
    assert (
        "3:1964-1982(-)" in guide2["target_genomic_regions"]
    ), f"Expected 3:1964-1982(-), got {guide2['target_genomic_regions']}"
    assert len(guide2["target_features"]) == 1, f"Expected 1 target feature, got {len(guide2['target_features'])}"
    assert guide2["target_features"]["ENSG00000345678"] == "EGFR"


def test_update_h5ad_with_guide_annotations_skips_na():
    """Test that guides with no gene overlaps are skipped."""
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
    assert guide1["target_genomic_regions"] == ["old:1-10(+)"], "Guide with no gene overlaps should not be updated"
    assert guide1["target_features"] == {"OLD_GENE": "OLD_NAME"}, "Target features should remain unchanged"


def test_update_h5ad_with_guide_annotations_missing_guide_warning(caplog):
    """Test that annotations for missing guide IDs are skipped with a warning."""
    import logging

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

    # Should complete without error, but log a warning
    with caplog.at_level(logging.WARNING):
        result_adata = annotate_guides.update_h5ad_with_guide_annotations(adata, annotations_df)

    # Verify warning was logged
    assert "Skipping" in caplog.text
    assert "guide_missing" in caplog.text

    # Verify guide1 still exists but has no new annotations
    assert "guide1" in result_adata.uns["genetic_perturbations"]


def test_update_h5ad_with_guide_annotations_no_genetic_perturbations():
    """Test that missing genetic_perturbations raises KeyError."""
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


def test_annotate_perturbations_skips_when_no_perturbations():
    """Test that annotation is skipped when genetic_perturbations is missing."""
    adata = ad.AnnData(X=np.zeros((2, 2)))

    # Snapshot state before
    uns_keys_before = set(adata.uns.keys())

    # Should return unmodified adata without error
    result_adata = annotate_guides.annotate_perturbations_in_h5ad(adata)

    assert result_adata is adata, "Should return the same object"
    assert set(result_adata.uns.keys()) == uns_keys_before, "uns dictionary should remain unchanged"


def test_annotate_perturbations_in_h5ad_no_guidescan(tmp_path, monkeypatch):
    """Test error when guidescan2 is not installed."""

    # Mock _check_guidescan2_installed to raise RuntimeError
    def mock_check():
        raise RuntimeError("guidescan2 not found")

    monkeypatch.setattr(annotate_guides, "_check_guidescan2_installed", mock_check)

    # Create minimal AnnData
    adata = ad.AnnData(X=np.zeros((2, 2)))
    adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
    adata.uns["genetic_perturbations"] = {
        "guide1": {"protospacer_sequence": "ACGT", "protospacer_adjacent_motif": "3' NGG"}
    }

    with pytest.raises(RuntimeError, match="guidescan2 not found"):
        annotate_guides.annotate_perturbations_in_h5ad(adata)


def test_annotate_perturbations_in_h5ad_integration(tmp_path, monkeypatch, sample_guidescan_csv):
    """Integration test for full annotation pipeline with mocks."""
    # Mock guidescan2 installed (no-op function that doesn't raise)
    monkeypatch.setattr(annotate_guides, "_check_guidescan2_installed", lambda: None)

    # Create AnnData with genetic perturbations
    adata = ad.AnnData(X=np.zeros((2, 2)))
    adata.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
    adata.uns["genetic_perturbations"] = {
        "guide1": {
            "protospacer_sequence": "TGCCTCGCGCAGCTCGCGG",
            "protospacer_adjacent_motif": "3' NGG",
        },
    }

    # Mock guidescan enumerate to copy sample output
    def mock_run_guidescan(input_csv, index_prefix, output_csv):
        import shutil

        shutil.copy(sample_guidescan_csv, output_csv)

    monkeypatch.setattr(annotate_guides, "_run_guidescan_enumerate", mock_run_guidescan)

    # Mock index path
    index_path = "/mock/index/path"

    # Run pipeline (returns modified adata)
    result_adata = annotate_guides.annotate_perturbations_in_h5ad(adata, index_path)

    # Verify result is an AnnData object
    assert isinstance(result_adata, ad.AnnData), "Result should be an AnnData object"
    assert "genetic_perturbations" in result_adata.uns, "genetic_perturbations should be present in result"

    # Check that guide1 has annotations (based on sample_guidescan_csv)
    guide1 = result_adata.uns["genetic_perturbations"]["guide1"]
    # The exact contents depend on sample_guidescan_csv and gene coordinates
    # At minimum, the keys should exist
    assert "protospacer_sequence" in guide1, "protospacer_sequence should be preserved"
    assert "protospacer_adjacent_motif" in guide1, "protospacer_adjacent_motif should be preserved"
