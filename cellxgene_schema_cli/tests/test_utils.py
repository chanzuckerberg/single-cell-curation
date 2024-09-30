import pandas as pd
import pytest
from anndata import AnnData
from cellxgene_schema.utils import (
    get_hash_digest_column,
    map_ontology_term,
    read_h5ad,
    remap_deprecated_features,
    remove_deprecated_features,
    replace_ontology_term,
)
from fixtures.examples_validate import adata, adata_non_raw, h5ad_valid


@pytest.fixture
def adata_with_raw():
    return adata.copy()


@pytest.fixture
def adata_without_raw():
    return adata_non_raw.copy()


@pytest.fixture
def deprecated_features():
    return ["ERCC-00002", "ENSG00000127603"]


@pytest.fixture
def remapped_features():
    return {"ENSMUSG00000059552": "ENSMUSG00000059552_NEW"}


@pytest.fixture
def deprecated_term_map_with_replacement_match():
    return {"EFO:0009899": "EFO:0000001"}


@pytest.fixture
def deprecated_term_map_no_replacement_match():
    return {"EFO:0000002": "EFO:0000003"}


def test_remove_deprecated_features__with_raw(adata_with_raw, deprecated_features):
    # Call the function under test
    result = remove_deprecated_features(adata=adata_with_raw, deprecated=deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]


def test_remove_deprecated_features__without_raw(adata_without_raw, deprecated_features):
    # Verify existing fixtures don't contain the deprecated features
    assert adata_without_raw.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552",
        "ENSSASG00005000004",
    ]

    # Call the function under test
    result = remove_deprecated_features(adata=adata_without_raw, deprecated=deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == ["ENSMUSG00000059552", "ENSSASG00005000004"]
    assert result.raw is None


def test_remap_deprecated_features__with_raw(adata_with_raw, remapped_features):
    # Verify existing fixtures don't contain the deprecated features
    assert adata_with_raw.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552",
        "ENSSASG00005000004",
    ]
    assert adata_with_raw.raw.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552",
        "ENSSASG00005000004",
    ]

    # Call the function under test
    result = remap_deprecated_features(adata=adata_with_raw, remapped_features=remapped_features)

    # Check if the deprecated features are replaced
    assert result.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552_NEW",
        "ENSSASG00005000004",
    ]
    assert result.raw.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552_NEW",
        "ENSSASG00005000004",
    ]


def test_remap_deprecated_features__without_raw(adata_without_raw, remapped_features):
    # Call the function under test
    result = remap_deprecated_features(adata=adata_without_raw, remapped_features=remapped_features)

    # Check if the deprecated features are replaced
    assert result.var_names.tolist() == [
        "ERCC-00002",
        "ENSG00000127603",
        "ENSMUSG00000059552_NEW",
        "ENSSASG00005000004",
    ]
    assert result.raw is None


def test_replace_ontology_term__with_replacement(adata_with_raw, deprecated_term_map_with_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_with_replacement_match)

    expected = ["EFO:0000001", "EFO:0008992"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
    assert sorted(actual) == expected


def test_replace_ontology_term__no_replacement(adata_with_raw, deprecated_term_map_no_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_no_replacement_match)
    expected = ["EFO:0008992", "EFO:0009899"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
    print(actual)
    assert all(a == b for a, b in zip(actual, expected))


def test_map_ontology_term__(adata_without_raw):
    update_map = {"donor_1": "CL:0000001", "donor_2": "CL:0000002"}
    map_ontology_term(adata_without_raw.obs, "cell_type", "donor_id", update_map)
    expected = ["CL:0000001", "CL:0000002"]
    actual = adata_without_raw.obs["cell_type_ontology_term_id"].dtype.categories
    assert all(a == b for a, b in zip(actual, expected))
    donor_1_rows = adata_without_raw.obs.loc[adata_without_raw.obs["donor_id"] == "donor_1"]
    assert all(a == "CL:0000001" for a in donor_1_rows["cell_type_ontology_term_id"])
    donor_2_rows = adata_without_raw.obs.loc[adata_without_raw.obs["donor_id"] == "donor_2"]
    assert all(a == "CL:0000002" for a in donor_2_rows["cell_type_ontology_term_id"])


class TestGetHashDigestColumn:
    def test_get_hash_digest_column(self, adata_with_raw):
        hash_digest_column = get_hash_digest_column(adata_with_raw.obs)
        expected_column = pd.Series(["ab6yl9v%fZ", "f-dZLjjiRl"], index=["X", "Y"])
        pd.testing.assert_series_equal(hash_digest_column, expected_column)


class TestReadH5AD:
    def test_read_h5ad(self):
        h5ad_path = h5ad_valid
        adata = read_h5ad(h5ad_path)
        assert isinstance(adata, AnnData)
        assert adata.isbacked
