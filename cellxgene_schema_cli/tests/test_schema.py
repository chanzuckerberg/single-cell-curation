import semver
from cellxgene_schema import schema


def test_get_current_schema_version():
    assert semver.Version.is_valid(schema.get_current_schema_version())
    assert schema.get_current_schema_version() == "5.2.0"
