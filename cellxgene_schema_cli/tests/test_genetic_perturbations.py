import tempfile

import pytest
from cellxgene_schema.validate import validate
from fixtures.examples_validate import (
    adata_gene_perturbations,
    adata_gene_perturbations_control,
    adata_gene_perturbations_invalid_bad_multi,
    adata_gene_perturbations_invalid_bad_strategy,
    adata_gene_perturbations_invalid_contains_derived,
    adata_gene_perturbations_invalid_contains_derived_features,
    adata_gene_perturbations_invalid_control_role_mismatch,
    adata_gene_perturbations_invalid_extra_key,
    adata_gene_perturbations_invalid_intended_features_bad_id,
    adata_gene_perturbations_invalid_intended_features_not_dict,
    adata_gene_perturbations_invalid_intended_features_versioned,
    adata_gene_perturbations_invalid_key_comma,
    adata_gene_perturbations_invalid_key_quote,
    adata_gene_perturbations_invalid_key_slash,
    adata_gene_perturbations_invalid_key_whitespace,
    adata_gene_perturbations_invalid_missing_key,
    adata_gene_perturbations_invalid_missing_obs_columns,
    adata_gene_perturbations_invalid_missing_pam,
    adata_gene_perturbations_invalid_missing_protospacer,
    adata_gene_perturbations_invalid_missing_role,
    adata_gene_perturbations_invalid_na_key,
    adata_gene_perturbations_invalid_obs_without_uns,
    adata_gene_perturbations_invalid_pam_format,
    adata_gene_perturbations_invalid_protospacer_chars,
    adata_gene_perturbations_invalid_protospacer_long,
    adata_gene_perturbations_invalid_protospacer_short,
    adata_gene_perturbations_invalid_role_enum,
    adata_gene_perturbations_invalid_strategy_without_id,
    adata_gene_perturbations_with_intended_features,
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
        adata_gene_perturbations_invalid_bad_strategy,
        adata_gene_perturbations_invalid_bad_multi,
        adata_gene_perturbations_invalid_missing_key,
        adata_gene_perturbations_invalid_control_role_mismatch,
        adata_gene_perturbations_invalid_contains_derived,
        adata_gene_perturbations_invalid_missing_obs_columns,
        adata_gene_perturbations_invalid_key_whitespace,
        adata_gene_perturbations_invalid_key_slash,
        adata_gene_perturbations_invalid_key_comma,
        adata_gene_perturbations_invalid_key_quote,
        adata_gene_perturbations_invalid_protospacer_short,
        adata_gene_perturbations_invalid_protospacer_long,
        adata_gene_perturbations_invalid_protospacer_chars,
        adata_gene_perturbations_invalid_pam_format,
        adata_gene_perturbations_invalid_missing_role,
        adata_gene_perturbations_invalid_missing_protospacer,
        adata_gene_perturbations_invalid_missing_pam,
        adata_gene_perturbations_invalid_role_enum,
        adata_gene_perturbations_invalid_contains_derived_features,
        adata_gene_perturbations_invalid_obs_without_uns,
        adata_gene_perturbations_invalid_na_key,
        adata_gene_perturbations_invalid_strategy_without_id,
        adata_gene_perturbations_invalid_extra_key,
        adata_gene_perturbations_invalid_intended_features_versioned,
        adata_gene_perturbations_invalid_intended_features_bad_id,
        adata_gene_perturbations_invalid_intended_features_not_dict,
    ],
)
def test_invalid_gene_perturbations(bad):
    success, errors = _validate_adata(bad)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"


def test_invalid_gene_perturbations_missing_obs_columns():
    """Test that missing obs columns when uns has genetic_perturbations is caught."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_missing_obs_columns)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    assert any(
        "When adata.uns['genetic_perturbations'] is present, obs must contain 'genetic_perturbation_id' and 'genetic_perturbation_strategy'."
        in error
        for error in errors
    )


def test_invalid_gene_perturbations_na_key():
    """uns['genetic_perturbations'] must not use 'na' as a perturbation identifier key."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_na_key)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    assert any("'na'" in error and "MUST NOT be present" in error for error in errors)


def test_invalid_gene_perturbations_strategy_without_id():
    """genetic_perturbation_strategy must not be present when genetic_perturbation_id is absent."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_strategy_without_id)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    assert any("genetic_perturbation_strategy" in error and "not present" in error for error in errors)


def test_invalid_gene_perturbations_extra_key():
    """Unknown extra keys in a genetic_perturbations entry must be rejected."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_extra_key)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    assert any("extra_annotation" in error and "MUST NOT be present" in error for error in errors)


def test_valid_gene_perturbations_with_intended_features():
    """intended_features with valid gene IDs (no version suffix) must pass validation."""
    success, errors = _validate_adata(adata_gene_perturbations_with_intended_features)
    assert success, errors


def test_invalid_intended_features_versioned_id():
    """intended_features key with ENS version suffix must fail (version must be stripped)."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_intended_features_versioned)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    # Versioned IDs are not found in the gene reference; the validator cannot infer the organism.
    assert any("ENSG00000141510.7" in error for error in errors)


def test_invalid_intended_features_bad_id():
    """intended_features key that is not a valid gene ID must fail validation."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_intended_features_bad_id)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
    assert any("NOT_A_GENE_ID" in error for error in errors)


def test_invalid_intended_features_not_dict():
    """intended_features value that is not a dict must fail validation."""
    success, errors = _validate_adata(adata_gene_perturbations_invalid_intended_features_not_dict)
    assert not success, f"Expected validation to fail but it succeeded. Errors: {errors}"
