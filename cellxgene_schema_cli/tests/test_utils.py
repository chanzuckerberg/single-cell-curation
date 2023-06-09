from unittest.mock import patch

import pytest
from cellxgene_schema.utils import get_deprecated_features


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
