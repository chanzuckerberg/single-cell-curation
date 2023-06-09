from unittest.mock import patch

import anndata as ad
import pytest
from cellxgene_schema.utils import get_deprecated_features, remove_deprecated_features
from fixtures.examples_validate import adata, adata_non_raw


@pytest.fixture
def organisms():
    return ["apple", "dog", "mouse"]


def test_get_deprecated_features(tmp_path, organisms):
    expected_deprecated_feature_ids = []
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            organism_feature_ids = [f"{organism}:{i}" for i in range(4)]
            for feature_id in organism_feature_ids:
                fp.write(feature_id + "\n")
            expected_deprecated_feature_ids.extend(organism_feature_ids)
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_features()
    expected_deprecated_feature_ids.sort()
    actual_deprecated_features.sort()
    assert expected_deprecated_feature_ids == actual_deprecated_features


def test_get_deprecated_features__no_files(tmp_path):
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_features()
    assert actual_deprecated_features == []


def test_get_deprecated_features__empty_feature_files(tmp_path, organisms):
    for organism in organisms:
        with open(f"{tmp_path}/{organism}_diff.txt", "w") as fp:
            fp.write("")
    with patch("cellxgene_schema.utils.env.ONTOLOGY_DIR", tmp_path):
        actual_deprecated_features = get_deprecated_features()
    assert actual_deprecated_features == []


def test_remove_deprecated_features():
    with patch("cellxgene_schema.utils.get_deprecated_features") as mock_get_deprecated_features:
        mock_get_deprecated_features.return_value = ["apple:0", "apple:1", "dog:0", "dog:1"]
        ad.AnnData()


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
    # Create a test AnnData object with some deprecated and non-deprecated features
    with patch("cellxgene_schema.utils.get_deprecated_features") as mock_get_deprecated_features:
        mock_get_deprecated_features.return_value = deprecated_features
        # Call the function under test
        result = remove_deprecated_features(adata_with_raw)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]


def test_remove_deprecated_features__without_raw(adata_without_raw, deprecated_features):
    # Create a test AnnData object with some deprecated and non-deprecated features
    with patch("cellxgene_schema.utils.get_deprecated_features") as mock_get_deprecated_features:
        mock_get_deprecated_features.return_value = deprecated_features
        # Call the function under test
        result = remove_deprecated_features(adata_without_raw)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw is None
