"""
Tests for guidescan2 gRNA annotation functionality.
"""

import shutil
import subprocess
from unittest.mock import MagicMock

import anndata as ad
import numpy as np
import pandas as pd
import pytest
from cellxgene_schema import annotate_guides


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
    input_csv.write_text(
        "id,sequence,pam,chromosome,start,end,sense\n"
        "test_guide,AAAAAAAAAAAAAAAAAAAAAGG,NGG,,,,\n"
    )

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
