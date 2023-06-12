import pytest
from cellxgene_schema.utils import remove_deprecated_features
from fixtures.examples_validate import adata, adata_non_raw


@pytest.fixture
def adata_with_raw():
    return adata.copy()


@pytest.fixture
def adata_without_raw():
    return adata_non_raw.copy()


@pytest.fixture
def deprecated_features():
    return ["ERCC-00002", "ENSG00000127603"]


def test_remove_deprecated_features__with_raw(adata_with_raw, deprecated_features):
    # Call the function under test
    result = remove_deprecated_features(adata_with_raw, deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]


def test_remove_deprecated_features__without_raw(adata_without_raw, deprecated_features):
    # Call the function under test
    result = remove_deprecated_features(adata_without_raw, deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw is None
