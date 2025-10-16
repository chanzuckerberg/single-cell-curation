"""Tests for HCA Validator."""

from pathlib import Path

import pytest
from hca_schema_validator import HCAValidator


def test_import():
    """Test that HCAValidator can be imported."""
    assert HCAValidator is not None


def test_instantiate():
    """Test that HCAValidator can be instantiated."""
    validator = HCAValidator()
    assert validator is not None
    assert hasattr(validator, 'validate_adata')


def test_inheritance():
    """Test that HCAValidator inherits from Validator."""
    from cellxgene_schema.validate import Validator
    
    validator = HCAValidator()
    assert isinstance(validator, Validator)


def test_has_validation_methods():
    """Test that HCAValidator has expected validation methods."""
    validator = HCAValidator()
    
    # Public methods
    assert hasattr(validator, 'validate_adata')
    
    # Protected methods (can override later)
    assert hasattr(validator, '_deep_check')
    assert hasattr(validator, '_validate_feature_ids')
    assert hasattr(validator, '_set_schema_def')


def test_validate_valid_h5ad():
    """Test validation of a valid h5ad file."""
    # Path to test fixture (relative to repo root)
    test_file = Path(__file__).parent.parent.parent / \
                "cellxgene_schema_cli/tests/fixtures/h5ads/example_valid.h5ad"
    
    if not test_file.exists():
        pytest.skip(f"Test fixture not found: {test_file}")
    
    validator = HCAValidator()
    is_valid = validator.validate_adata(str(test_file))
    
    # Should pass validation
    assert is_valid is True
    assert len(validator.errors) == 0
    # May have warnings (that's okay)
    assert validator.warnings is not None


def test_validate_invalid_h5ad():
    """Test validation of an invalid h5ad file."""
    # Path to test fixture (relative to repo root)
    test_file = Path(__file__).parent.parent.parent / \
                "cellxgene_schema_cli/tests/fixtures/h5ads/example_invalid_CL.h5ad"
    
    if not test_file.exists():
        pytest.skip(f"Test fixture not found: {test_file}")
    
    validator = HCAValidator()
    is_valid = validator.validate_adata(str(test_file))
    
    # Should fail validation
    assert is_valid is False
    assert len(validator.errors) > 0
