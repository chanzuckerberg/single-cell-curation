from unittest.mock import patch

import semver
from cellxgene_schema import schema


# Mock the version to ensure consistent test results
@patch("cellxgene_schema.schema.__version__", "6.0.0")
def test_get_current_schema_version():
    assert semver.Version.is_valid(schema.get_current_schema_version())
    assert schema.get_current_schema_version() == "6.0.0"
