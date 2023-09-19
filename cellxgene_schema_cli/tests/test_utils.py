import numpy as np
import pytest
from anndata import AnnData
from cellxgene_schema.utils import (
    enforce_canonical_format,
    map_ontology_term,
    read_h5ad,
    remove_deprecated_features,
    replace_ontology_term,
)
from fixtures.examples_validate import adata, adata_non_raw, h5ad_valid
from scipy.sparse import coo_matrix


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
def deprecated_term_map_with_replacement_match():
    return {"EFO:0009899": "EFO:0000001"}


@pytest.fixture
def deprecated_term_map_no_replacement_match():
    return {"EFO:0000002": "EFO:0000003"}


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


def test_replace_ontology_term__with_replacement(adata_with_raw, deprecated_term_map_with_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_with_replacement_match)

    expected = ["EFO:0009918", "EFO:0000001"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
    assert all(a == b for a, b in zip(actual, expected))


def test_replace_ontology_term__no_replacement(adata_with_raw, deprecated_term_map_no_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_no_replacement_match)
    expected = ["EFO:0009899", "EFO:0009918"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
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


@pytest.fixture
def noncanonical_matrix():
    array = np.array([[1, 0, 1], [3, 2, 3], [4, 5, 4]])
    return coo_matrix((array[0], (array[1], array[2])))


class TestEnforceCanonical:
    def test_adata_with_noncanonical_X_and_raw_X(self, noncanonical_matrix):
        assert noncanonical_matrix.has_canonical_format is False
        adata = AnnData(noncanonical_matrix)
        enforce_canonical_format(adata)
        assert adata.X.has_canonical_format is True

    def test_adata_with_noncanonical_raw_X(self, noncanonical_matrix):
        assert noncanonical_matrix.has_canonical_format is False
        adata = AnnData(raw=AnnData(noncanonical_matrix))
        enforce_canonical_format(adata)
        assert adata.raw.X.has_canonical_format is True

    def test_adata_with_canonical_raw_X(self, adata_with_raw):
        enforce_canonical_format(adata)
        assert adata_with_raw.raw.X.has_canonical_format is True

    def test_adata_with_canonical_X(self, adata_without_raw):
        enforce_canonical_format(adata)
        assert adata_without_raw.X.has_canonical_format is True


class TestReadH5AD:
    def test_read_h5ad(self):
        h5ad_path = h5ad_valid
        adata = read_h5ad(h5ad_path)
        assert isinstance(adata, AnnData)
        assert adata.isbacked

    def test_read_h5ad_to_memory(self):
        # Provide a valid h5ad path or a valid object resembling a path
        h5ad_path = h5ad_valid
        adata = read_h5ad(h5ad_path, to_memory=True)
        assert isinstance(adata, AnnData)
        assert not adata.isbacked
