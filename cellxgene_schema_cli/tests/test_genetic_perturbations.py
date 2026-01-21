import tempfile

import pytest
from cellxgene_schema.validate import validate
from fixtures.examples_validate import (
    adata_gene_perturbations,
    adata_gene_perturbations_control,
    adata_gene_perturbations_invalid_bad_multi,
    adata_gene_perturbations_invalid_bad_strategy,
    adata_gene_perturbations_invalid_contains_derived,
    adata_gene_perturbations_invalid_contains_na,
    adata_gene_perturbations_invalid_control_role_mismatch,
    adata_gene_perturbations_invalid_missing_key,
    adata_gene_perturbations_invalid_missing_obs_columns,
)


def _validate_adata(adata):
    with tempfile.TemporaryDirectory() as tmp:
        p = tmp + "/temp.h5ad"
        adata.copy().write_h5ad(p, compression="gzip")
        success, errors, _ = validate(p)
        return success, errors


def test_valid_gene_perturbations():
    success, errors = _validate_adata(adata_gene_perturbations)
    assert success, errors


def test_valid_gene_perturbations_control():
    success, errors = _validate_adata(adata_gene_perturbations_control)
    assert success, errors


@pytest.mark.parametrize(
    "bad",
    [
        adata_gene_perturbations_invalid_contains_na,
        adata_gene_perturbations_invalid_bad_strategy,
        adata_gene_perturbations_invalid_bad_multi,
        adata_gene_perturbations_invalid_missing_key,
        adata_gene_perturbations_invalid_control_role_mismatch,
        adata_gene_perturbations_invalid_contains_derived,
        adata_gene_perturbations_invalid_missing_obs_columns,
    ],
)
def test_invalid_gene_perturbations(bad):
    success, errors = _validate_adata(bad)
    assert not success


def test_invalid_gene_perturbations_missing_obs_columns():
    """Test that missing obs columns when uns has genetic_perturbations is caught."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_missing_obs_columns)
    assert not success
    assert any(
        "When adata.uns['genetic_perturbations'] is present, obs must contain 'genetic_perturbation_id' and 'genetic_perturbation_strategy'."
        in error
        for error in errors
    )
