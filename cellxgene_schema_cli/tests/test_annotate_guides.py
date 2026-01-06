"""
Tests for guidescan2 gRNA annotation functionality.
"""

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
