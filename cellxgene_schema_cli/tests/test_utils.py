import pandas as pd
import pytest
from anndata import AnnData
from cellxgene_schema.utils import (
    get_chunks,
    get_hash_digest_column,
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
