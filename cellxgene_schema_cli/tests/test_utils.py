import numpy as np
import pandas as pd
import pytest
from anndata import AnnData
from cellxgene_schema.utils import (
    get_chunks,
    get_hash_digest_column,
    getattr_anndata,
    map_ontology_term,
    move_column_from_obs_to_uns,
    read_h5ad,
    remap_deprecated_features,
    remove_deprecated_features,
    replace_delimiter,
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
    return ["ENSG00000141510", "ENSG00000127603"]


@pytest.fixture
def remapped_features():
    return {"ENSG00000012048": "ENSG00000012048_NEW"}


@pytest.fixture
def deprecated_term_map_with_replacement_match():
    return {"EFO:0009899": "EFO:0000001"}


@pytest.fixture
def deprecated_term_map_no_replacement_match():
    return {"EFO:0000002": "EFO:0000003"}


def test_remove_deprecated_features__with_raw(adata_with_raw, deprecated_features):
    # Verify existing fixtures don't contain the deprecated features
    assert adata_with_raw.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]

    # Call the function under test
    result = remove_deprecated_features(adata=adata_with_raw, deprecated=deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == [
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]
    assert result.raw.var_names.tolist() == [
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]


def test_remove_deprecated_features__without_raw(adata_without_raw, deprecated_features):
    # Verify existing fixtures don't contain the deprecated features
    assert adata_without_raw.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]

    # Call the function under test
    result = remove_deprecated_features(adata=adata_without_raw, deprecated=deprecated_features)

    # Check if the deprecated features are removed
    assert result.var_names.tolist() == [
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]
    assert result.raw is None


def test_remap_deprecated_features__with_raw(adata_with_raw, remapped_features):
    # Verify existing fixtures don't contain the deprecated features
    assert adata_with_raw.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]
    assert adata_with_raw.raw.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]

    # Call the function under test
    result = remap_deprecated_features(adata=adata_with_raw, remapped_features=remapped_features)

    # Check if the deprecated features are replaced
    assert result.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048_NEW",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]
    assert result.raw.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048_NEW",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]


def test_remap_deprecated_features__without_raw(adata_without_raw, remapped_features):
    # Call the function under test
    result = remap_deprecated_features(adata=adata_without_raw, remapped_features=remapped_features)

    # Check if the deprecated features are replaced
    assert result.var_names.tolist() == [
        "ENSG00000127603",
        "ENSG00000141510",
        "ENSG00000012048_NEW",
        "ENSG00000139618",
        "ENSG00000002330",
        "ENSG00000000005",
        "ENSG00000000419",
    ]
    assert result.raw is None


def test_replace_ontology_term__with_replacement(adata_with_raw, deprecated_term_map_with_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_with_replacement_match)

    expected = ["EFO:0000001"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
    assert sorted(actual) == expected


def test_replace_ontology_term__no_replacement(adata_with_raw, deprecated_term_map_no_replacement_match):
    replace_ontology_term(adata_with_raw.obs, "assay", deprecated_term_map_no_replacement_match)
    expected = ["EFO:0009899"]
    actual = adata_with_raw.obs["assay_ontology_term_id"].dtype.categories
    print(actual)
    assert all(a == b for a, b in zip(actual, expected))


def test_map_ontology_term(adata_without_raw):
    update_map = {"donor_1": "CL:0000001", "donor_2": "CL:0000002"}
    map_ontology_term(adata_without_raw.obs, "cell_type", "donor_id", update_map)
    expected = ["CL:0000001", "CL:0000002"]
    actual = adata_without_raw.obs["cell_type_ontology_term_id"].dtype.categories
    assert all(a == b for a, b in zip(actual, expected))
    donor_1_rows = adata_without_raw.obs.loc[adata_without_raw.obs["donor_id"] == "donor_1"]
    assert all(a == "CL:0000001" for a in donor_1_rows["cell_type_ontology_term_id"])
    donor_2_rows = adata_without_raw.obs.loc[adata_without_raw.obs["donor_id"] == "donor_2"]
    assert all(a == "CL:0000002" for a in donor_2_rows["cell_type_ontology_term_id"])


def test_move_column_from_obs_to_uns(adata_with_raw):
    assert "assay_ontology_term_id" in adata_with_raw.obs.columns
    assert "assay_ontology_term_id" not in adata_with_raw.uns

    move_column_from_obs_to_uns(adata_with_raw, "assay_ontology_term_id")

    assert "assay_ontology_term_id" not in adata_with_raw.obs.columns
    assert adata_with_raw.uns["assay_ontology_term_id"] == "EFO:0009899"


def test_replace_delimiter(adata_with_raw):
    adata_with_raw.obs["self_reported_ethnicity_ontology_term_id"] = "HsapDv:0000003,HsapDv:0000004"

    replace_delimiter(adata_with_raw.obs, ",", " || ", "self_reported_ethnicity_ontology_term_id")

    assert adata_with_raw.obs["self_reported_ethnicity_ontology_term_id"].eq("HsapDv:0000003 || HsapDv:0000004").all()


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


class TestChunks:
    def test_get_chunks__even(self):
        chunks = get_chunks(step_size=100, total_size=2000)
        assert len(chunks) == 20
        assert chunks[0] == (0, 100)
        assert chunks[-1] == (1900, 2000)
        for chunk_start, chunk_end in chunks:
            assert chunk_start < chunk_end
            assert chunk_end - chunk_start == 100

    def test_get_chunks__uneven(self):
        chunks = get_chunks(step_size=100, total_size=1850)
        assert len(chunks) == 19
        assert chunks[0] == (0, 100)
        assert chunks[-1] == (1800, 1850)

    def test_get_chunks__slightly_higher(self):
        chunks = get_chunks(step_size=100, total_size=2001)
        assert len(chunks) == 21
        assert chunks[0] == (0, 100)
        assert chunks[-1] == (2000, 2001)

    def test_get_chunks__total_less_than_step(self):
        chunks = get_chunks(step_size=100, total_size=20)
        assert len(chunks) == 1
        assert chunks[0] == (0, 20)


class TestGetattrAnndata:
    """Test suite for getattr_anndata function with nested access support."""

    @pytest.fixture
    def adata_with_nested_data(self):
        """Create an AnnData object with nested structures for testing."""
        obs = pd.DataFrame(
            {
                "cell_type": ["T cell", "B cell", "NK cell"],
                "assay_ontology_term_id": ["EFO:0009899", "EFO:0009899", "EFO:0009899"],
            },
            index=["cell1", "cell2", "cell3"],
        )
        var = pd.DataFrame(
            {
                "gene_symbol": ["GENE1", "GENE2", "GENE3"],
                "feature_type": ["Gene Expression", "Gene Expression", "Gene Expression"],
            },
            index=["ENSG000001", "ENSG000002", "ENSG000003"],
        )
        X = np.zeros((obs.shape[0], var.shape[0]), dtype=np.float32)
        uns = {
            "species_ontology_term_id": "NCBITaxon:9606",
            "some_list": [{"key1": "value1", "key2": "value2"}, {"key1": "value3", "key2": "value4"}],
            "nested_dict": {"level1": {"level2": "deep_value"}},
            "simple_list": ["a", "b", "c"],
        }
        adata = AnnData(X=X, obs=obs, var=var, uns=uns)
        # Add raw data (reuse same data structures)
        adata.raw = AnnData(X=X.copy(), obs=obs.copy(), var=var.copy())
        return adata

    def test_simple_attribute_access(self, adata_with_nested_data):
        """Test accessing simple attributes like obs, var, uns."""
        assert getattr_anndata(adata_with_nested_data, "obs") is adata_with_nested_data.obs
        assert getattr_anndata(adata_with_nested_data, "var") is adata_with_nested_data.var
        assert getattr_anndata(adata_with_nested_data, "uns") is adata_with_nested_data.uns
        assert getattr_anndata(adata_with_nested_data, "X") is adata_with_nested_data.X

    @pytest.mark.parametrize(
        "component,column,expected_column",
        [
            ("obs", "cell_type", "cell_type"),
            ("obs", "assay_ontology_term_id", "assay_ontology_term_id"),
            ("var", "gene_symbol", "gene_symbol"),
        ],
    )
    def test_dataframe_column_access(self, adata_with_nested_data, component, column, expected_column):
        """Test accessing columns in obs and var."""
        result = getattr_anndata(adata_with_nested_data, f"{component}.{column}")
        expected = getattr(adata_with_nested_data, component)[expected_column]
        pd.testing.assert_series_equal(result, expected)

    def test_empty_string_attribute(self, adata_with_nested_data):
        """Test that empty string returns None."""
        result = getattr_anndata(adata_with_nested_data, "")
        assert result is None

    def test_multiple_consecutive_dots(self, adata_with_nested_data):
        """Test that multiple consecutive dots raise an exception."""
        with pytest.raises(ValueError, match="Invalid attribute path"):
            getattr_anndata(adata_with_nested_data, "obs..cell_type")

    def test_uns_dict_access(self, adata_with_nested_data):
        """Test accessing dictionary keys in uns."""
        result = getattr_anndata(adata_with_nested_data, "uns.species_ontology_term_id")
        assert result == "NCBITaxon:9606"

    @pytest.mark.parametrize(
        "index,expected",
        [
            (0, "a"),
            (1, "b"),
            (2, "c"),
        ],
    )
    def test_uns_list_indexing(self, adata_with_nested_data, index, expected):
        """Test accessing list elements in uns using bracket notation."""
        result = getattr_anndata(adata_with_nested_data, f"uns.simple_list[{index}]")
        assert result == expected

    @pytest.mark.parametrize(
        "index,expected",
        [
            (0, {"key1": "value1", "key2": "value2"}),
            (1, {"key1": "value3", "key2": "value4"}),
        ],
    )
    def test_uns_nested_list_dict_access(self, adata_with_nested_data, index, expected):
        """Test accessing nested dictionaries within lists."""
        result = getattr_anndata(adata_with_nested_data, f"uns.some_list[{index}]")
        assert result == expected

    def test_uns_nested_dict_access(self, adata_with_nested_data):
        """Test accessing nested dictionaries."""
        result = getattr_anndata(adata_with_nested_data, "uns.nested_dict")
        assert result == {"level1": {"level2": "deep_value"}}

    def test_raw_access(self, adata_with_nested_data):
        """Test accessing raw data."""
        result = getattr_anndata(adata_with_nested_data, "raw.var")
        assert result is adata_with_nested_data.raw.var

        # raw doesn't have obs attribute, should return None
        result = getattr_anndata(adata_with_nested_data, "raw.obs")
        assert result is None

    def test_raw_obs_column_access(self, adata_with_nested_data):
        """Test that raw.obs doesn't exist (raw only has var, not obs)."""
        # raw.obs should return None since raw doesn't have obs
        result = getattr_anndata(adata_with_nested_data, "raw.obs.cell_type")
        assert result is None

    def test_none_attribute(self, adata_with_nested_data):
        """Test that None attribute returns None."""
        result = getattr_anndata(adata_with_nested_data, None)
        assert result is None

    def test_missing_key_in_uns(self, adata_with_nested_data):
        """Test accessing non-existent key in uns returns None."""
        result = getattr_anndata(adata_with_nested_data, "uns.non_existent_key")
        assert result is None

    def test_missing_column_in_obs(self, adata_with_nested_data):
        """Test accessing non-existent column in obs returns None."""
        result = getattr_anndata(adata_with_nested_data, "obs.non_existent_column")
        assert result is None

    def test_list_index_out_of_bounds(self, adata_with_nested_data):
        """Test accessing list index out of bounds returns None."""
        result = getattr_anndata(adata_with_nested_data, "uns.simple_list[10]")
        assert result is None

    def test_list_index_negative(self, adata_with_nested_data):
        """Test that negative indices are not supported (returns None)."""
        result = getattr_anndata(adata_with_nested_data, "uns.simple_list[-1]")
        assert result is None

    def test_adata_without_raw(self):
        """Test accessing raw when it doesn't exist returns None."""
        obs = pd.DataFrame({"cell_type": ["T cell"]}, index=["cell1"])
        var = pd.DataFrame({"gene": ["GENE1"]}, index=["ENSG000001"])
        X = np.zeros((obs.shape[0], var.shape[0]), dtype=np.float32)
        adata = AnnData(X=X, obs=obs, var=var)
        adata.raw = None

        result = getattr_anndata(adata, "raw.var")
        assert result is None

    @pytest.mark.parametrize(
        "path,expected",
        [
            ("uns.some_list[0].key1", "value1"),
            ("uns.some_list[1].key2", "value4"),
        ],
    )
    def test_complex_nested_path(self, adata_with_nested_data, path, expected):
        """Test complex nested paths."""
        result = getattr_anndata(adata_with_nested_data, path)
        assert result == expected

    def test_multiple_list_indices(self, adata_with_nested_data):
        """Test accessing nested lists."""
        adata_with_nested_data.uns["nested_list"] = [[1, 2, 3], [4, 5, 6]]
        first_list = getattr_anndata(adata_with_nested_data, "uns.nested_list[0]")
        assert first_list == [1, 2, 3]
        assert first_list[1] == 2
