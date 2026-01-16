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
    ],
)
def test_invalid_gene_perturbations(bad):
    success, errors = _validate_adata(bad)
    assert not success
